#!/usr/bin/env python
from __future__ import division
from past.utils import old_div
import os
import glob
import EMAN2
import datetime
import logger
import mpi
import random
import argparse
import shutil
import EMAN2db
import numpy as np
import tqdm
import sys
import subprocess
import ctypes
try:
    from sphire.libpy import sp_global_def
    from sphire.libpy import sp_utilities
    from sphire.libpy import sp_fundamentals
    from sphire.libpy import sp_applications
    from sphire.libpy import sp_morphology
    from sphire.libpy import sp_filter
    from sphire.libpy import sp_projection
except ImportError:
    import sp_global_def
    import sp_utilities
    import sp_fundamentals
    import sp_applications
    import sp_morphology
    import sp_filter
    import sp_projection
    
# Set default values (global variables written in ALL CAPS)
TIMESTAMP_LENGTH= 23  # chars
FILTRAD_DEFAULT= 0.5  # px^-1
MAPTHRESH_SIGMA= 1.0  # sigma
MASKFILT_SPFREQ= .02  # px^-1
MASKTHRESH_DEF= 0.05  # 0..1

# AVGFILT mode
VOL_AVG_UNFILT=      'vol_avg_unfil.hdf'  # Average of input half-set maps, unfiltered
VOL_AVG_FILT=        'vol_avg_fil.hdf'    # Low-pass filtered average of input half-set maps

# SP_MASK mode
VOL_FILTERED=        'input_filtered.hdf'  # Low-pass filtered input map
MASK_FILE=           'mask_final.hdf'      # Final mask
VOL_MASKED=          'vol_masked.hdf'      # Diagnostic masked input map

# PROJSUBTRACT mode
CTF_PARAMS=          'ctf_params.txt'      # CTF parameter file
MPI_STACK_DIR=       'Stacks'              # Directory for MPI stacks
MPI_PROJ_PREFIX=     'proj'                # Stem of reprojection stacks, BDBs numbered, MRCSs preceded by 'stk', MPI threads have process number
COMB_PROJ_BDB_NAME=  'all_projections'     # Combined BDB of projections, which will be subtracted from input images
DOC_DIR=             'Docs'                # Directory for text files
SELECT_DOC=          'select'              # Stem of selection doc file, numbered with MPI-thread process number
MPI_ORIG_PREFIX=     'orig'                # Stem of copies of original stacks, numbered with MPI-thread process number
MPI_SUBTRACT_PREFIX= 'subtract'            # Stem of signal-subtracted stacks, BDBs numbered, MRCSs preceded by 'stk', MPI threads have process number
NORM_STATS_DOC=      'docnormstats.txt'    # Normalization statistics
COMB_SUBTRACT_BDB=   'all_subtracted'      # Combined BDB stack of signal-subtracted images
TEST_MONTAGE=        'stktestmontage.hdf'  # Montage of test images

# CENTERSHIFT mode
INITIAL_PARAM_DOC=   'params004_to_be_centered.txt'  # Original parameters
CENTERED_PARAM_DOC=  'params005_centered.txt'        # Parameters after centering

USAGE= """ 
PURPOSE:
Performs various steps involved in signal subtraction.

General syntax:
%s subtraction_mode --verbosity=2
Parameters:
  subtraction_mode : avgfilt, sp_mask, projsubtract, or centershift
  --verbosity : Verbosity level (0..3) (default: 2)

To sum and filter two maps:
%s avgfilt vol_0_unfil_027.hdf --avol2 vol_1_unfil_027.hdf --outdir Outdir --apix 1.34 --filtrad 13.4 
Parameters:
  --avol1 : Input map #1 (required)
  --avol2 : Input map #2 (optional)
  --outdir : Output directory
  --apix : Pixel size, Angstroms/pixel
  --filtrad : Low-pass filter radius, Angstroms (or if, apix not provided, pixels^-1)
Outputs:
  vol_avg_unfil.hdf : Unfiltered average map
  vol_avg_fil.hdf : Low-pass filtered average map

To generate an adaptive mask:
%s sp_mask --map2mask map_to_make_mask_of.mrc --fullmap vol_avg_unfil.hdf --outdir Outdir --mapthresh 0.005 --maskthresh 0.02 --masklowpass 13.4 --apix 1.34 --falloff 0.01
Inputs:
  --map2mask : Map from which to make a mask
  --fullmap : Map to multiply by mask
Parameters:
  --mapthresh : Threshold at which the initial map is binarized (default: +1 s.d.)
  --ndilation : Number of dilations (default: 3)
  --edge_width : Edge falloff, volxels (default: 5)
  --masklowpass : Radius of low-pass filter applied to input map, Angstroms (or if, apix not provided, pixels^-1)
  --falloff : Low-pass filter falloff in absolute frequency (default: 0.01 px^-1)
  --allow_disconnected : Allow disconnected regions in the mask (default: False)
Outputs:
  mask_final.hdf : Final, adaptive mask
  masked_vol.hdf : Masked input map

To generate projections of a reconstruction and subtract from experimental images:
%s projsubtract --origparts bdb:Particles#data --map2subtract map_to_subtract.hdf --projparams final_params_027.txt --outdir Outdir --normalize --nmontage=12 --inmem --saveprojs --stats
Inputs:
  --origparts : Input original particles
  --map2subtract : Map whose projections are to be subtracted
  --projparams : Alignment parameters
Outputs:
  ctf_params.txt : CTF parameters for images in input stack
  bdb:Outdir#to_be_subtracted : Re-projections of map to be subtracted
  bdb:Outdir#subtracted : Subtracted stack
Parameters:
  --normalize : Normalizes images (recommended)
Advanced parameters:
  --nmontage : Saves examples of original/projected/subtracted images
  --inmem : Stores projections in memory (slightly faster, but far more memory-intensive)
  --saveprojs : Saves re-projections
  --stats : Writes average and sigma for each particle before normalization

To center reconstruction:
%s centershift --cvol1 vol_0_unfil_027.hdf --cvol2 vol_1_unfil_027.hdf --shiftparams final_params_027.txt --diffimgs bdb:Outdir#subtracted --outdir Outdir --volradius 52 --prefix='vol'
Inputs:
  --cvol1 : Input map #1
  --cvol2 : Input map #2 (optional)
  --shiftparams : Alignment parameters for each image (optional, if not specified will read from header)
  --diffimgs : Difference images (optional, shifted parameters will be written here)
Parameters:
  --apix : Pixel size, Angstroms/pixel
  --volradius : Radius (Angstroms or, if apix not provided, pixels) of the centered blob to which your reconstruction will be aligned
  --prefix : Prefix for output reconstructions (default: vol)
Outputs:
  vol001_average.hdf : Average of inputs cvol1 and cvol2, i.e., reconstruction of signal-subtracted images
  vol002_inverted.hdf : Negative of above reconstruction, subtracted region shouldn't have significant densities
  vol003_init_centered.hdf : Centered reconstruction
  params004_to_be_centered.txt : Copy of original alignment parameters
  params005_centered.txt : Alignment parameters after centering

""" % ((__file__,)*5)

MODIFIED="Modified 2021-01-26"

"""
CHANGELOG:
    2021-01-09 (trs) -- normalization during subtraction on by default
    2021-01-05 (trs & ab) -- intermediate re-projections are deleted by default
    2021-01-05 (trs) -- comparison images can be written
    2021-01-05 (trs & ab) -- to save memory, doesn't store entire output stack in memory
    2020-12-09 (trs) -- updated for Python 3
    2020-06-04 (trs) -- implemented MPI for subtraction step
    2020-03-05 (trs & ab) -- implemented MPI for reprojection step
    2019-09-27 (trs) -- mask can be disconnected
    2019-07-17 (trs) -- reconstruction outsourced to MERIDIEN
    2019-05-31 (trs) -- can low-pass filter map before masking
    2019-05-22 (trs) -- replaced soft mask with sp_mask
    2019-05-22 (trs) -- updated for SPHIRE v1.2
    TODO: don't use safe-exit() with errors
    TODO: add full set of options for sp_mask
"""

class BdbNames:
    """
    Adapted from sp_signalsubtract.py
    
    Class containing the following attributes:
        bdb_stem : Stem, no directory, no leading "bdb:", no trailing ".bdb"
        bdb_dir : Directory, not counting "EMAN2DB"
        em2db_dir : Directory, including "EMAN2DB"
        bdb_name : Of the form "bdb:DirectoryName#stem"
        bdb_path : Of the form "DirectoryName/EMAN2DB/stem.bdb"
    """

    def __init__(self, bdb_stem, dir_above_eman2db=None):
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
        
def setup_mpi(verbose=False):
    """ 
    Sets up MPI.
    
    Argument:
        verbose : (boolean) whether to print to screen
    Returns:
        mpidict : (dictionary) {'use', 'myid', 'main_process', 'is_main', 'size'}
    """
    
    main_process= 0  # default
    
    usempi= "OMPI_COMM_WORLD_SIZE" in os.environ
    
    if usempi:
        number_of_proc= mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        myid= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        
        #  choose a random node as main
        if myid == 0: main_process= random.randint(0,number_of_proc-1)
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
            log= logger.Logger(
                base_logger=logger.BaseLogger_Files(), 
                base_logger2=logger.BaseLogger_Print(), 
                file_name=logname
                )
            verbose= False  # logger output will be echoed to screen
        else:
            log= logger.Logger(base_logger=logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        if is_main: print("WARNING: Using old logger.py library")
        log= logger.Logger(base_logger=logger.BaseLogger_Files())#, file_name=logname)
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
        mesg : message to write
        log : instance of Logger class
        verbose : (boolean) whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    if is_main:
        if verbose: print(mesg)
        if log: log.add(mesg)

def quick_barrier():
    """
    Synchronizes parallel processes before continuing
    Safe for non-MPI calls also
    """
    
    if "OMPI_COMM_WORLD_SIZE" in os.environ : mpi.mpi_barrier(mpi.MPI_COMM_WORLD) 	
    
def safe_exit(mesg=None, verbosity=0):
    """
    Properly closes MPI before exiting, preventing ugly warnings
    Safe for non-MPI calls also
    
    Modified 2021-01-17

    Argument:
        verbosity : controls how much information to write to screen (0..3)
    """
    
    if "OMPI_COMM_WORLD_SIZE" in os.environ: 
        my_rank= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        if verbosity >= 3:
            print("Process %s: Synchronizing with barrier" % my_rank)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD) 
        if verbosity >= 2:
            print("Process %s: Barrier reached" % my_rank)
        mpi.mpi_finalize()
    else:
        my_rank= 0
    
    if verbosity >= 1:
        print("Process %s: Finalized" % my_rank)
    
    if mesg:
        sp_global_def.ERROR(mesg, __file__, 1)
    else:
        sp_global_def.print_timestamp( "Finish" )
    
    exit()
    
def printvars(variables, typeTF=False, quitTF=False):
    """
    Print the local variables in the caller's frame.
    
    Adapted from https://stackoverflow.com/questions/6618795/get-locals-from-calling-namespace-in-python
    """
    
    import inspect
    import six

    if type(variables) is list:
        # Weird things happen if 
        assert isinstance(variables[0], six.string_types), "UH OH! Passed non-string %s instead of variable name" % variables[0]
        
        variable_list= variables
    elif isinstance(variables, six.string_types):  # six works with both Pythons 2 & 3
        variable_list= [variables]
    else:
        print("ERROR!! Don't know how to deal with type %s" % type(variables) )
        exit()
    
    frame = inspect.currentframe()
    dictionary= frame.f_back.f_locals
    
    print("")
    for variable in variable_list : 
        try: 
            mesg= "%s: %s" % (variable, dictionary[variable])
            if typeTF : mesg+= " %s" % type(dictionary[variable])
            print(mesg)
        except KeyError: 
            print('ERROR!! printvars')
            print(dictionary)
            print('EXITING!!')
            exit()
        
    del frame
    
    if quitTF: 
        print('\nExiting printvars...')  # reminder in case I forget to take out the quit flag
        exit()

def setup_subprocess():
    """
    Determines Python version
    Subprocess usage differs for Python versions 2 & 3, as does the syntax for /dev/null
    
    Returns
        Python version number
        allowed instance of DEVNULL (specific to Python version)
        
    Copied from sp_eval_isac.py, 2020-06-01
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
    
def signalsubtract(options, outdir='.', verbosity=0):
    """
    Main function overseeing various signal-subtraction modes.
    
    Arguments:
        options : (namepspace) command-line options, run 'sp_signalsubtract.py -h' for an exhaustive list
        outdir : output directory
        verbosity : controls how much information to write to screen (0..3)
    """

    # Set up MPI
    mpioptions= setup_mpi()
    is_main=mpioptions['is_main']
    
    # Set output directory and log file name
    log_obj, verbose= prepare_outdir_log(outdir, verbosity>=1, is_main=is_main)
    quick_barrier()
    sp_global_def.write_command(outdir)  # assumes mpi_comm_rank=0, which isn't necessarily main process, which is random
    quick_barrier()
    
    # Determine mode
    if options.mode == 'avgfilt':
        # If no volume provided, give error
        if not options.avol1 and not options.avol2:
            safe_exit("No volumes provided!")
        else:
            # If only one volume provided, use that one
            if bool(options.avol1) ^ bool(options.avol2):  # "^" is XOR
                if options.avol2: avol1= avol2= options.avol2
                if options.avol1: avol2= avol1= options.avol1
            # Both volumes provided
            else:
                avol1= options.avol1
                avol2= options.avol2
        
        if options.filtrad:
            # If pixel size provided, convert to reciprocal pixels
            if options.apix:
                filtrad= old_div(options.apix, options.filtrad)
            else:
                filtrad= options.filtrad
            
            if verbosity>=1 : print_log_msg(f"Will filter to {filtrad: .3f} px^-1", log_obj)
        else:
            # If no filter radius specified, simply filter to Nyquist
            filtrad= FILTRAD_DEFAULT
            if verbosity>=1 : print_log_msg(f"No filtrad specified, will filter to {filtrad: .3f} px^-1", log_obj)
        
        avgfilt(
            avol1, 
            avol2, 
            os.path.join(outdir, VOL_AVG_UNFILT), 
            os.path.join(outdir, VOL_AVG_FILT), 
            filtrad, 
            log=log_obj, 
            verbose=verbosity>=1, 
            is_main=mpioptions['is_main']
            )

    elif options.mode == 'sp_mask':
        continueTF= True  # Will proceed unless some information is missing
        
        if not options.fullmap:
            if verbosity>=1 : print_log_msg("No fullmap given, will not multiply by mask", log_obj)
            # UNTESTED
        
        if not options.map2mask:
            continueTF= False
            safe_exit("No map2mask given!")
        
        if continueTF:
            map2mask= options.map2mask
            
            if options.masklowpass:
                filt_falloff= options.falloff
                
                # If pixel size provided, convert to reciprocal pixels
                if options.apix:
                    filter_radius_px= old_div(options.apix, options.masklowpass)
                else:
                    filter_radius_px= options.masklowpass
                
                mesg= "Will filter to %s px^-1 with a falloff of %s" % (filter_radius_px, filt_falloff)
                if verbosity>=1 : print_log_msg(mesg, log_obj)
            
            if options.mapthresh:
                mapthresh= options.mapthresh
                if verbosity>=1 : print_log_msg("Will use a mapthresh of %s" % mapthresh, log_obj)
            else:
                vol= sp_utilities.get_image(map2mask)
                mapthresh= vol['mean'] + vol['sigma']
                mesg= "No mapthresh given, will use %s (mean + 1*sigma)" % mapthresh
                if verbosity>=1 : print_log_msg(mesg, log_obj)
            
        if continueTF:
            spmask(
                options.fullmap, 
                map2mask, 
                mapthresh, 
                options.ndilation, 
                options.edge_width, 
                os.path.join(outdir, MASK_FILE), 
                os.path.join(outdir, VOL_MASKED), 
                filtfile=os.path.join(outdir, VOL_FILTERED), 
                filtrad=filter_radius_px, 
                falloff=filt_falloff, 
                log=log_obj, 
                verbose=verbosity>=1, 
                is_main=mpioptions['is_main']
                )
        else:
            exit()

    elif options.mode == 'projsubtract':
        continueTF= True  # Will proceed unless some information is missing
        
        # Filenames/parameters which should have been supplied
        if options.origparts:
            if is_main and verbosity>=1 : 
                mesg= "Will subtract projections from stack '%s'" % options.origparts
                print_log_msg(mesg, log_obj)
        else:
            continueTF= False
            safe_exit("No origparts given!")
            
        if options.map2subtract:
            if is_main and verbosity>=1 : 
                mesg= "Will store projections from map '%s'" % options.map2subtract
                print_log_msg(mesg, log_obj)
        else:
            continueTF= False
            safe_exit("No map2subtract given!")
            
        if options.projparams:
            if is_main and verbosity>=1 : 
                mesg= "Will use projection parameters from file '%s'" % options.projparams
                print_log_msg(mesg, log_obj)
        else:
            continueTF= False
            safe_exit("No projparams given!")
            
        if options.nonorm:
            if is_main and verbosity>=1 : 
                mesg= "Will skip normalization of input images before subtraction"
                print_log_msg(mesg, log_obj)
        else:
            options.stats= True
            
            if is_main and verbosity>=1 : 
                mesg= "Will normalize input images before subtraction"
                print_log_msg(mesg, log_obj)
            
        if options.normalize:
            if options.nonorm : 
                continueTF= False
                safe_exit("Options '--normalize' and '--nonorm' are incompatible")

            if is_main and verbosity>=1 : 
                mesg= "DEPRECATION WARNING: Option '--normalize no longer needed', is now default"
                print_log_msg(mesg, log_obj)
        
        if options.stats: 
            stats_doc= os.path.join(outdir, NORM_STATS_DOC)
            if is_main and verbosity>=1: 
                mesg= "Will write input image statistics"
                print_log_msg(mesg, log_obj)
        else:
            stats_doc= None
                
        if is_main and verbosity>=1 : 
            if options.nmontage!=0:
                mesg= "Will write montage of up to %s test images" % options.nmontage
                print_log_msg(mesg, log_obj)
        
            # Memory estimate
            num_imgs= EMAN2.EMUtil.get_image_count(options.origparts)
            img= sp_utilities.get_im(options.origparts, 0)
            mem_gb= num_imgs*img['nx']*img['ny']*4.0/1000000000.0
            mesg= f"Will require{mem_gb: .3f} GB memory" 
            
            if options.inmem:
                mesg+= " (total), will compute re-projections in memory"
            else:
                mesg+= ", will write re-projections incrementally (see '--inmem')"
            print_log_msg(mesg, log_obj)
        
            if options.saveprojs:
                mesg= "Will keep re-projections\n"
            else:
                mesg= "Will delete re-projections after subtraction\n"
            print_log_msg(mesg, log_obj)
        
        if continueTF : 
            projsubtract(
                options.origparts, 
                os.path.join(outdir, CTF_PARAMS), 
                options.map2subtract, 
                options.projparams, 
                outdir, 
                mpioptions, 
                reprojs=COMB_PROJ_BDB_NAME, 
                diffimg_stem=COMB_SUBTRACT_BDB, 
                stats_doc=stats_doc, 
                do_norm=(not options.nonorm), 
                numtestimgs=options.nmontage,
                inmem=options.inmem, 
                saveprojs=options.saveprojs, 
                log=log_obj, 
                verbosity=verbosity, 
                )
        else:
            exit()

    elif options.mode == 'centershift':
        continueTF= True  # Will proceed unless some information is missing
        
        # If no volume provided, give error
        if not options.cvol1 and not options.cvol2:
            safe_exit("No volumes provided!")
            continueTF= False
        else:
            # If only one volume provided, use that one
            if bool(options.cvol1) ^ bool(options.cvol2):  # "^" is XOR
                if options.cvol2: cvol1= cvol2= options.cvol2
                if options.cvol1: cvol2= cvol1= options.cvol1
            # Both volumes provided
            else:
                cvol1= options.cvol1
                cvol2= options.cvol2
        
        if options.shiftparams:
            inparams= options.shiftparams
            if verbosity>=1:
                mesg= "Will apply alignment parameters to '%s'" % options.shiftparams
                print_log_msg(mesg, log_obj)
        else:
            inparams= None
        
        if options.diffimgs:
            mesg= "No alignment parameters given, will extract from %s" % options.diffimgs
            diffimgs= options.diffimgs
        else:
            diffimgs= None
            if not inparams:
                if verbosity>=1: 
                    mesg= "No difference-image stack nor alignment parameters given, not writing parameters file"
                    print_log_msg(mesg, log_obj)
        
        if not options.volradius:
            continueTF= False
            safe_exit("No volradius given!")
        else:
            # Convert to pixels if pixel size is provided
            if options.apix:
                volradiuspx= old_div(options.volradius, options.apix)
                if verbosity>=1:
                    mesg= "Will use a radius of %spx for centering (%sA at %s A/px)" % (volradiuspx, options.volradius, options.apix)
                    print_log_msg(mesg, log_obj)
            else:
                volradiuspx= options.volradius
                if verbosity>=1:
                    mesg= "No pixel size provided, will use a radius of %spx for centering" % volradiuspx
                    print_log_msg(mesg, log_obj)
            
        if continueTF:
            centershift(
                cvol1, 
                cvol2, 
                diffimgs, 
                os.path.join(outdir, options.prefix + '001_average.hdf'), 
                os.path.join(outdir, options.prefix + '002_inverted.hdf'), 
                volradiuspx, 
                os.path.join(outdir, options.prefix + '003_centered.hdf'), 
                inparams, 
                os.path.join(outdir, INITIAL_PARAM_DOC), 
                os.path.join(outdir, CENTERED_PARAM_DOC), 
                log_obj, 
                verbosity>=1, 
                is_main
                )
        else:
            exit()
        
    else:
        safe_exit("No valid mode given! Acceptable modes are: avgfilt, sp_mask, projsubtract, centershift")

    if is_main : print_log_msg("Done!", log_obj)
    
def avgfilt(input1, input2, outunfilt, outfilt, filtrad=0.5, log=None, verbose=False, is_main=True):
    """
    Averages two maps and applies low-pass filter.
    
    Arguments:
        input1 : input map #1
        input2 : input map #2
        outunfilt : output averaged map filename
        outfilt : output low-pass filtered map filename
        filtrad : (float) low-pass filter radius, reciprocal pixels
        log : instance of Logger class
        verbose : (boolean) -- whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    # Average
    outvol= average2(input1, input2, outunfilt, log, verbose, is_main)
    
    # Low-pass filter
    if is_main and verbose:
        print_log_msg(f"Filtering '{outunfilt}' to {filtrad: .3f}px^-1", log)
        cmd= "e2proc3d.py %s --process=filter.lowpass.gauss:cutoff_freq=%s %s" % (outunfilt, filtrad, outfilt)
        print_log_msg(cmd, log)
    
    outvol.process_inplace('filter.lowpass.gauss', {'cutoff_freq':filtrad})
    outvol.write_image(outfilt)
    
    if is_main and verbose:
        print_log_msg(f"Wrote filtered map to '{outfilt}'\n", log)

def average2(input1, input2, outfile, log, verbose, is_main):
    """
    Averages two images or volumes.
    
    Arguments:
        input1 : input #1
        input2 : input #2
        outfile : output average filename
        log : instance of Logger class
        verbose : (boolean) -- whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    Returns:
        outavg : EMData object of average
    """
    
    cmd= "e2proc3d.py --addfile %s %s --mult=0.5 %s" % (input1, input2, outfile)
    if verbose and is_main:
        print_log_msg("Averaging '%s' and '%s'" % (input1, input2), log)
        print_log_msg(cmd, log)
    outavg= sp_utilities.get_image(input1)
    outavg.add( sp_utilities.get_image(input2) )
    outavg.mult(0.5)
    outavg.write_image(outfile)
    if verbose and is_main:
        print_log_msg("Wrote averaged map to '%s'\n" % outfile, log)
    
    return outavg

def spmask(fullmapfile, premaskfile, mapthresh, ndilation, edgewidth, maskfile, maskedfile, 
            filtfile=None, filtrad=None, falloff=0.01, log=None, verbose=False, is_main=True):
    """
    Generates soft mask and applies to input volume.
    
    Arguments:
        fullmapfile : input map to be multiplied by mask (optional)
        premaskfile : input map from which the mask will be generated
        mapthresh : (float) map threshold for initial binarization
        ndilation : (int) number of dilations
        edgewidth : (int) falloff of smooth edge, voxels
        maskfile : output final mask
        maskedfile : output map after multiplication by final mask
        filtfile : optional low-pass filtered pre-mask
        filtrad : (float) low-pass filter radius, absolute frequency
        falloff : (float) filter falloff, absolute frequency
        log : instance of Logger class
        verbose : (boolean) -- whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    # Make mask
    if is_main:
        if verbose:
            print_log_msg(
                "Thresholding input map '%s' at a density of %s with %s dilations and edge width %s to '%s'" % 
                (premaskfile, mapthresh, ndilation, edgewidth, maskfile), log)
            cmd= "sp_mask.py %s %s --threshold='%s' --ndilation=%s --edge_width=%s\n" % (premaskfile, maskfile, mapthresh, ndilation, edgewidth)
            print_log_msg(cmd, log)
        
        # Read volume as EMData object
        input_vol= sp_utilities.get_image(premaskfile)
        
        # Optional low-pass filter
        if filtrad:
            input_vol= sp_filter.filt_tanl(input_vol, filtrad, falloff)
            input_vol.write_image(filtfile)
        
        maskvol= sp_morphology.adaptive_mask_scipy(
            input_vol,
            threshold=mapthresh,
            ndilation=ndilation,
            edge_width=edgewidth,
            do_print=True,
            )
        
        # Write mask to disk
        maskvol.write_image(maskfile)
        if verbose : print_log_msg("Wrote mask to '%s'" % maskfile, log)
        
        # Multiply by mask
        fullmapobj= sp_utilities.get_image(fullmapfile)
        if verbose:
            mesg= "Multiplying map '%s' by mask '%s'" % (fullmapfile, maskfile)
            print_log_msg(mesg, log)
            cmd= "e2proc3d.py %s %s --multfile=%s" % (fullmapfile, maskedfile, maskfile)
            print_log_msg(cmd, log)
        fullmapobj.mult(maskvol)
        fullmapobj.write_image(maskedfile)
        if verbose : print_log_msg("Wrote mask-multiplied map to '%s'\n" % maskedfile, log)
    
def projsubtract(parts0, ctfdoc, map2subtract, projparams, outdir, mpioptions, 
                reprojs=COMB_PROJ_BDB_NAME, diffimg_stem=COMB_SUBTRACT_BDB, stats_doc=None, 
                do_norm=True, numtestimgs=0, inmem=False, saveprojs=False, 
                log=None, verbosity=0):
    """
    Computes re-projections of map, and subtracts them from image stack.
    
    Arguments:
        parts0 : input starting image stack from which re-projections will be subtracted
        ctfdoc : output CTF parameter doc
        map2subtract : input map from which re-projections will be computed
        projparams : input alignment parameter doc
        outdir : output directory
        reprojs : output re-projection stackname (no "bdb:stackname" or "stackname.bdb", simply "stackname" )
        diffimg_stem : output difference stackname (no "bdb:stackname" or "stackname.bdb", simply "stackname" )
        stats_doc : optional output file with average and sigma of input images
        do_norm : (boolean) perform normalization of input images
        numtestimgs : number of montages of test comparisons to write
        inmem : (boolean) whether to retain re-projections in RAM
        saveprojs : (boolean) whether to save or delete re-projections
        log : instance of Logger class
        verbosity : how much information to write to screen (0..3)
        mpioptions : (dict) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    # BDB names
    reproj_bdb_class= BdbNames(reprojs, outdir)
    diff_bdb_class= BdbNames(diffimg_stem, outdir)
    diffimg_stack= os.path.join(outdir, diffimg_stem + '.mrcs')
    
    # Filenames for parallelization
    mpi_stack_dir= os.path.join(outdir, MPI_STACK_DIR)
    mpi_orig_template= 'bdb:' + outdir + '#' + MPI_ORIG_PREFIX + '_mpi{:03d}'
    mpi_proj_bdb_stem= MPI_PROJ_PREFIX + '_mpi{:03d}'
    mpi_partdoc_template= os.path.join(outdir, DOC_DIR, SELECT_DOC + '_mpi{:03d}.txt')
    mpi_diff_bdb_stem= MPI_SUBTRACT_PREFIX + '_mpi{:03d}'
    
    if mpioptions['use']:
        my_mpi_proc_id= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        mpi_mrc_stack_template= os.path.join(mpi_stack_dir, 'stk' + mpi_proj_bdb_stem + '.mrcs') 
        mpi_mrc_stack_file= mpi_mrc_stack_template.format(my_mpi_proc_id)
    else:
        my_mpi_proc_id= 0
        mpi_mrc_stack_template= os.path.join(mpi_stack_dir, 'stk' + MPI_PROJ_PREFIX + '.mrcs') 
        mpi_mrc_stack_file= mpi_mrc_stack_template
    is_main= mpioptions['is_main']  # not as useful as my_mpi_proc_id
    
    if is_main:
        # Export CTF parameters and compute re-projections of map to be subtracted
        if verbosity>=1 : 
            print_log_msg("Exporting CTF parameters...", log)
            cmd= "sp_header.py %s --params=ctf --export=%s" % (parts0, ctfdoc)
            print_log_msg(cmd, log)
        sp_applications.header(parts0, 'ctf', fexport=ctfdoc)
        if verbosity>=1 : 
            print_log_msg("Wrote CTF parameters to '%s'\n" % ctfdoc, log)
    
        # Re-project map to be subtracted
        if verbosity>=1 : 
            mesg= "Projecting map '%s' using angles in '%s' and CTF parameters '%s'" %(
                map2subtract, projparams, ctfdoc)
            print_log_msg(mesg, log)
            cmd= "sp_project3d.py %s %s --angles=%s --CTF=%s --tril\n" % (map2subtract, reproj_bdb_class.bdb_name, projparams, ctfdoc)
            print_log_msg(cmd, log)
    
    quick_barrier()
    
    # Compute projections
    project3d_MPI(
        map2subtract, 
        projparams, 
        outdir, 
        ctfdoc, 
        parts0, 
        mpi_partdoc_template, 
        reproj_bdb_class, 
        mpioptions, 
        mpi_orig_template=mpi_orig_template, 
        mpi_proj_bdb_stem=mpi_proj_bdb_stem, 
        mpi_stack_dir=mpi_stack_dir, 
        mpi_mrc_stack_file=mpi_mrc_stack_file, 
        do_inmem=inmem, 
        log=log, 
        verbosity=verbosity
        )
    
    quick_barrier()
    
    # Subtract re-projections from images
    if is_main and verbosity>=1:
        mesg= "Subtracting re-projections '%s' from original images '%s'" % (reproj_bdb_class.bdb_name, parts0)
        print_log_msg(mesg, log)
        cmd= "sp_process.py %s %s %s --subtract_stack" % (parts0, reproj_bdb_class.bdb_name, diff_bdb_class.bdb_name)
        if do_norm : cmd+= ' --normalize'
        print_log_msg(cmd, log)
    
    quick_barrier()
    
    subtract_stack_MPI(
        parts0, 
        reproj_bdb_class.bdb_name, 
        diffimg_stem, 
        outdir, 
        stats_doc=stats_doc, 
        do_norm=do_norm, 
        mpi_orig_template=mpi_orig_template, 
        mpi_proj_bdb_stem=mpi_proj_bdb_stem, 
        mpi_stack_dir=mpi_stack_dir, 
        mpi_diff_bdb_stem=mpi_diff_bdb_stem, 
        mpioptions=mpioptions, 
        log=log, 
        verbosity=verbosity
        )
    
    quick_barrier()
    
    if is_main:
        if numtestimgs!=0:
            if verbosity>=1:
                print_log_msg("Writing up to %s example images" % numtestimgs, log)
            
            ####write_example_montages(
                ####mpi_orig_template, 
                ####mpi_proj_bdb_stem, 
                ####mpi_diff_bdb_stem, 
                ####numtestimgs,
                ####outdir, 
                ####TEST_MONTAGE,
                ####log=log, 
                ####verbosity=verbosity
                ####)
            write_example_montages(
                parts0, 
                reproj_bdb_class.bdb_name, 
                diff_bdb_class.bdb_name, 
                numtestimgs,
                outdir, 
                TEST_MONTAGE,
                log=log, 
                verbosity=verbosity
                )
            
        if not saveprojs:
            if verbosity>=1:
                print_log_msg("Cleaning up intermediate projections", log)
            
            clean_projs(
                reproj_bdb_class.bdb_name,
                mpi_proj_bdb_stem, 
                mpi_mrc_stack_template, 
                outdir, 
                mpioptions, 
                inmem=inmem, 
                log=log, 
                verbosity=verbosity
                )
        else:
            print_log_msg("Keeping projections", log)
    
    quick_barrier()

def project3d_MPI(volfile, projparams, outdir, ctfdoc, orig_part_bdb, 
                  mpi_partdoc_template, comb_proj_bdb_obj, mpioptions, 
                  mpi_orig_template=None, mpi_proj_bdb_stem=None, 
                  mpi_stack_dir=None, mpi_mrc_stack_file=None, 
                  do_inmem=False, log=None, verbosity=0):
    """
    Computes a series of projections.
    
    Arguments:
        volfile : 3D volume from which projections will be generated
        projparams : input projection parameters file
        outdir : output directory
        ctfdoc : input CTF parameter file
        comb_proj_bdb_obj : output instance of class BdbNames for projection stack
        mpi_partdoc_template : template for output particle-selection files
        mpioptions : (dict) MPI options
        orig_part_bdb : BDB stack to be split up for MPI
        mpi_orig_template : template for virtual stacks of original BDB split up for MPI
        mpi_proj_bdb_stem : BDB stem (w/o leading 'bdb:') for output projection stacks
        mpi_stack_dir : directory where output MRC stacks will be written
        mpi_mrc_stack_file : MRC stack file (if inmem)
        do_inmem : (boolean) compute projection in RAM
        log : instance of Logger class
        verbosity : how much information to write to screen (0..3)
    """
    
    # Check if using MPI
    RUNNING_UNDER_MPI= "OMPI_COMM_WORLD_SIZE" in os.environ
    main_mpi_proc= 0
    is_main= mpioptions['is_main']
    
    # TODO: Validate that files needed for MPI are defined/present
    
    # Set up MPI
    if RUNNING_UNDER_MPI:
        my_mpi_proc_id= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        num_mpi_procs   = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        sp_global_def.MPI= True
        
        # If running on more than one node, one process per node will load shared data
        node_communicator= mpi.mpi_comm_split_type(mpi.MPI_COMM_WORLD, mpi.MPI_COMM_TYPE_SHARED,  0, mpi.MPI_INFO_NULL)
        my_mpi_proc_on_node_id= mpi.mpi_comm_rank(node_communicator)
        
        python_version, devnull= setup_subprocess()
        if my_mpi_proc_id == main_mpi_proc:
            where_stdout= None
        else:
            where_stdout= devnull
    else:
        my_mpi_proc_id= 0
        num_mpi_procs= 1
        my_mpi_proc_on_node_id= 0
        node_communicator= 0
        
    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()
        if verbosity>=1 and my_mpi_proc_id==main_mpi_proc : print_log_msg("Disabling BDB cache", log)
    
    # Make directories (only for main process)
    if my_mpi_proc_id == main_mpi_proc:
        # Set up output directories & files
        if not os.path.isdir(outdir) : os.makedirs(outdir)
        bdb_dir= os.path.join(outdir, 'EMAN2DB')
        if not os.path.isdir(bdb_dir) : os.makedirs(bdb_dir)
        if not os.path.isdir(mpi_stack_dir) : os.makedirs(mpi_stack_dir)
        doc_dir= os.path.dirname(mpi_partdoc_template)
        if not os.path.isdir(doc_dir) : os.makedirs(doc_dir)
            
    angles_list= sp_utilities.read_text_row(projparams)
    ctf_list= sp_utilities.read_text_row(ctfdoc)
    start, end= sp_applications.MPI_start_end( len(angles_list), num_mpi_procs, my_mpi_proc_id )
    num_angles= end - start + 1
    
    # Read specific attribute as list, rather than read image and then get attribute
    coord_list= EMAN2.EMUtil.get_all_attributes(orig_part_bdb, "ptcl_source_coord")
    source_list= EMAN2.EMUtil.get_all_attributes(orig_part_bdb, "ptcl_source_image")
    # (Reading the image individual and then the attributes took ~30 minutes in test set.)
    
    quick_barrier()
    
    if RUNNING_UNDER_MPI:
        disp_unit= np.dtype("f4").itemsize
        
        # Each main process on each main node will 
        if my_mpi_proc_on_node_id == 0:
            if my_mpi_proc_id==main_mpi_proc and verbosity>=1:
                print_log_msg("Preparing map '%s' for projection" % volfile, log)
            invol= sp_utilities.get_im(volfile)
            idim= invol['ny']
            volft= sp_projection.prep_vol(invol, npad= 2, interpolation_method= 1)
            if my_mpi_proc_id==main_mpi_proc and verbosity>=1:
                print_log_msg('Finished preparing map, sharing memory for parallelization...\n', log)
            numpy_volft= EMAN2.EMNumPy.em2numpy(volft)
            nxvol= volft.get_xsize()
            nyvol= volft.get_ysize()
            nzvol= volft.get_zsize()
            size= nxvol * nyvol * nzvol
        else:
            nxvol= 0
            nyvol= 0
            nzvol= 0
            size= 0
            idim= 0
            
        # Broadcast dimensions to all processes
        nxvol= sp_utilities.bcast_number_to_all(nxvol, source_node= main_mpi_proc, mpi_comm= node_communicator)
        nyvol= sp_utilities.bcast_number_to_all(nyvol, source_node= main_mpi_proc, mpi_comm= node_communicator)
        nzvol= sp_utilities.bcast_number_to_all(nzvol, source_node= main_mpi_proc, mpi_comm= node_communicator)
        idim = sp_utilities.bcast_number_to_all(idim,  source_node= main_mpi_proc, mpi_comm= node_communicator)
        
        # Allocate memory
        win_sm, base_ptr= mpi.mpi_win_allocate_shared(size*disp_unit, disp_unit, mpi.MPI_INFO_NULL, node_communicator)
        
        if my_mpi_proc_on_node_id != main_mpi_proc:
            base_ptr,= mpi.mpi_win_shared_query(win_sm, mpi.MPI_PROC_NULL)
            
        size= nxvol * nyvol * nzvol
        """
        The following no longer works:
        volbuffer= np.frombuffer(np.core.multiarray.int_asbuffer(base_ptr, size*disp_unit), dtype= 'f4')
        Using Adnan's workaround instead
        """
        ptr_n= ctypes.cast(base_ptr, ctypes.POINTER( ctypes.c_int * size) )
        volbuffer= np.frombuffer(ptr_n.contents, dtype="f4")
        
        volbuffer= volbuffer.reshape(nzvol, nyvol, nxvol)

        if( my_mpi_proc_on_node_id == main_mpi_proc ):
            np.copyto(volbuffer, numpy_volft)
            del numpy_volft, invol, volft

        emnumpy2= EMAN2.EMNumPy()
        volft= emnumpy2.register_numpy_to_emdata(volbuffer)
        volft.set_attr_dict({'is_complex':1,  'is_complex_x': 0, 'is_fftodd': 0, 'is_fftpad': 1, 'is_shuffled': 1,'npad': 2})
    
        # Generate list of indices (e2bdb will need it)
        index_list= list( range(start, end) )  # will work in both Python 2 & 3
        selection_file= mpi_partdoc_template.format(my_mpi_proc_id)
        sp_utilities.write_text_row(index_list, selection_file)
        
        ##### Read specific attribute as list, rather than read image and then get attribute
        ####coord_list= EMAN2.EMUtil.get_all_attributes(orig_part_bdb, "ptcl_source_coord")
        ####source_list= EMAN2.EMUtil.get_all_attributes(orig_part_bdb, "ptcl_source_image")
        ##### (Reading the image individual and then the attributes took ~30 minutes in test set.)
    
        # Set up BDB
        mpi_proj_bdb_obj= BdbNames(mpi_proj_bdb_stem.format(my_mpi_proc_id), outdir)
        mpi_bdb_name= mpi_proj_bdb_obj.bdb_name
        mpi_bdb_path= mpi_proj_bdb_obj.bdb_path
        
        # Only for one process
        if start == 0:
            print_log_msg("Projecting map '%s' in %s globs of about %s images" % (volfile, num_mpi_procs, end), log)
    
    # Not using MPI
    else:
        invol= sp_utilities.get_im(volfile)
        nxvol= invol['nx']
        nyvol= idim= invol['ny']
        print_log_msg("Preparing volume '%s' for projection" % volfile, log)
        volft= sp_projection.prep_vol(invol, npad= 2, interpolation_method= 1)
        print_log_msg('Finished preparing volume\n', log)
        
        # Set up BDB
        mpi_bdb_path= os.path.join(comb_proj_bdb_obj.bdb_path)
        mpi_bdb_name= comb_proj_bdb_obj.bdb_name
    
    # Initialize
    if os.path.exists(mpi_bdb_path) : os.remove(mpi_bdb_path)  # will otherwise merge with pre-existing file
    quick_barrier()
    part_counter= 0
    
    if do_inmem:
        new_bdb_dict= EMAN2db.db_open_dict(mpi_bdb_name)
        aligned_stack_obj= EMAN2.EMData(idim, idim, num_angles-1)
    
    # Show progress bar only for one process (or none)
    disableTF= start!=0 or verbosity<=1
    
    # Loop through angles
    for global_num in tqdm.tqdm(range(start, end), unit='proj', disable=disableTF, file=sys.stdout):
        angles= angles_list[global_num]
        prj_obj= sp_projection.prgl(volft, [angles[0], angles[1], angles[2], -angles[3], -angles[4] ], interpolation_method=1, return_real=False)
        sp_utilities.set_params_proj(prj_obj, angles[0:5])
        
        # Apply CTF
        ctf= sp_utilities.generate_ctf([
            ctf_list[global_num][0], 
            ctf_list[global_num][1], 
            ctf_list[global_num][2], 
            ctf_list[global_num][3], 
            ctf_list[global_num][4], 
            ctf_list[global_num][5], 
            ctf_list[global_num][6], 
            ctf_list[global_num][7]
            ])
        prj_obj.set_attr_dict({"is_complex": 0})  # should already be 0
        EMAN2.Util.mulclreal(prj_obj, sp_morphology.ctf_img_real(prj_obj.get_ysize(), ctf))
        prj_obj.set_attr_dict({"padffted":1, "is_complex":1})
        prj_obj= sp_fundamentals.fft(prj_obj)
        prj_obj.set_attr("ctf", ctf)
        prj_obj.set_attr("ctf_applied", 0)
        
        # Get image attributes
        if do_inmem:
            prj_dict= prj_obj.get_attr_dict()
            prj_dict["ptcl_source_coord_id"]= part_counter
            prj_dict["data_path"]= sp_utilities.makerelpath(comb_proj_bdb_obj.eman2db_dir, mpi_mrc_stack_file)
            
            # For verification
            ####if RUNNING_UNDER_MPI:
            prj_dict["ptcl_source_coord"]= coord_list[global_num]
            prj_dict["ptcl_source_image"]= source_list[global_num]
            
            new_bdb_dict[part_counter]= prj_dict
            aligned_stack_obj.insert_clip(prj_obj, (0, 0, part_counter) )
        else:
            prj_obj["ptcl_source_coord_id"]= part_counter
        
            # For verification
            ####if RUNNING_UNDER_MPI:
            prj_obj["ptcl_source_coord"]= coord_list[global_num]
            prj_obj["ptcl_source_image"]= source_list[global_num]
        
            prj_obj.write_image(mpi_bdb_name, part_counter)
        
        part_counter+= 1
    # End angles loop
    
    if do_inmem: aligned_stack_obj.write_image(mpi_mrc_stack_file)
    EMAN2db.db_close_dict(mpi_bdb_name)
    
    # Split input stack so that subtraction can be parallelized
    if RUNNING_UNDER_MPI:
        mpi_orig_name= mpi_orig_template.format(my_mpi_proc_id)
        cmd=  "e2bdb.py"
        args= "%s --makevstack %s --list %s" % (orig_part_bdb, mpi_orig_name, selection_file)
        if start==0 and verbosity>=1:
            print_log_msg("Finished making projections\n", log)
            print_log_msg("Making substacks from input '%s' for each parallel process" % orig_part_bdb, log)
            print_log_msg(cmd + " " + args, log)
        
        system_call_23(cmd, args, stdout=where_stdout)
        
        if start==0 and verbosity>=1:
            print_log_msg("Waiting for parallel processes to finish...\n", log)
    
    # Synchronize
    quick_barrier()
    
    # Write combined BDB (main process only)
    if RUNNING_UNDER_MPI and my_mpi_proc_id == main_mpi_proc:
        if verbosity>=1 : print_log_msg("Finished waiting, now combining stack files", log)
        
        # Combine stacks
        merge_bdbs(
            outdir, 
            mpi_proj_bdb_stem, 
            num_mpi_procs, 
            comb_proj_bdb_obj, 
            verbosity=verbosity, 
            log=log
            )
        
    quick_barrier()
    
def subtract_stack_MPI(minuend_stack, subtrahend_stack, combined_diffimg_bdb_stem, outdir, 
                       stats_doc=None, do_norm=False, mpi_orig_template=None, 
                       mpi_proj_bdb_stem=None, mpi_stack_dir=None, mpi_diff_bdb_stem=None, 
                       mpioptions=None, log=None, verbosity=0):
    """
    Subtracts one image stack from another, adapted from sp_process.py, 2018-12-06.
    
    Arguments:
        minuend_stack : input image stack from which subtrahend images will be subtracted
        subtrahend_stack : input image stack to be subtracted
        combined_diffimg_bdb_stem : BDB stem (w/o leading 'bdb:') of output difference-image stack
        outdir : output directory
        stats_doc : output text file with average and sigma of original images before subtraction
        do_norm : (boolean) perform normalization
        mpi_orig_template : template for virtual stacks of original BDB split up for MPI
        mpi_proj_bdb_stem : BDB stem (w/o leading 'bdb:') for input projection stacks
        mpi_stack_dir : directory where output MRC stacks will be written
        mpi_diff_bdb_stem : BDB stem (w/o leading 'bdb:') for output difference-image stacks
        mpioptions : (dict) MPI options
        log : instance of Logger class
        verbosity : how much information to write to screen (0..3)
    """
    
    # Check if using MPI
    RUNNING_UNDER_MPI= "OMPI_COMM_WORLD_SIZE" in os.environ
    main_mpi_proc= mpioptions['main_process']
    is_main= mpioptions['is_main']
    num_init_parts= EMAN2.EMUtil.get_image_count(minuend_stack)
    
    # Set up MPI
    if RUNNING_UNDER_MPI:
        my_mpi_proc_id= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        num_mpi_procs   = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        mpi_minuend_stack= mpi_orig_template.format(my_mpi_proc_id)
        mpi_subtrahend_stack= 'bdb:' + outdir + '#' + mpi_proj_bdb_stem.format(my_mpi_proc_id)
        mpi_difference_stack= os.path.join(mpi_stack_dir, 'stk' + mpi_diff_bdb_stem.format(my_mpi_proc_id) + '.mrcs') 
        mpi_diff_bdb_obj= BdbNames(mpi_diff_bdb_stem.format(my_mpi_proc_id), outdir)
    else:
        my_mpi_proc_id= 0
        num_mpi_procs= 1
        mpi_minuend_stack= minuend_stack
        mpi_subtrahend_stack= subtrahend_stack
        mpi_difference_stack= os.path.join(mpi_stack_dir, 'stk' + MPI_SUBTRACT_PREFIX + '.mrcs') 
        mpi_diff_bdb_obj= BdbNames(COMB_SUBTRACT_BDB, outdir)
        
    # Sanity checks
    try:
        nimages= EMAN2.EMUtil.get_image_count(mpi_minuend_stack)
    except: 
        printvars(['mpi_minuend_stack', 'my_mpi_proc_id'])
        safe_exit("Can't process mpi_minuend_stack")
    try:
        mimages= EMAN2.EMUtil.get_image_count(mpi_subtrahend_stack)
    except: 
        printvars(['mpi_subtrahend_stack', 'my_mpi_proc_id'])
        safe_exit("Can't process mpi_subtrahend_stack")
    assert nimages == mimages, 'ERROR!! Two input stacks have different number of images %s %s' % (nimages, mimages)
    
    new_bdb_dict= EMAN2db.db_open_dict(mpi_diff_bdb_obj.bdb_name)
    
    disableTF= my_mpi_proc_id !=0 or verbosity<=1
    
    for idx in tqdm.tqdm(range(nimages), unit='img', disable=disableTF, file=sys.stdout):
        try:
            mimage = sp_utilities.get_im(mpi_minuend_stack, idx)
        except: 
            printvars(['mpi_minuend_stack', 'idx'])
            safe_exit("Can't process mpi_minuend_stack")
            # (NOTE: I once got an error here which I wasn't able to reproduce.)
        
        try:
            simage= sp_utilities.get_im(mpi_subtrahend_stack, idx)
        except:
            printvars(['mpi_subtrahend_stack', 'idx'])
            safe_exit("Can't process mpi_subtrahend_stack")
            
        # Record statistics
        if stats_doc or do_norm:
            if idx == 0:
                radius= mimage.get_xsize()//2 - 1 
                mask= sp_utilities.model_circle(radius, mimage.get_xsize(), mimage.get_ysize())
            st = EMAN2.Util.infomask(mimage, mask, False)
        
        # Normalize
        if do_norm:
            mimage-= st[0]
            mimage/= st[1]
        
        diffimg= EMAN2.Util.subn_img(mimage, simage)
        ctf_obj= mimage.get_attr('ctf')
        diffimg.set_attr('ctf_applied', 0)
        diffimg.set_attr('ctf', ctf_obj)
        
        prj_obj= simage.get_attr('xform.projection')
        diffimg.set_attr('xform.projection', prj_obj)
        diffimg.write_image(mpi_difference_stack, idx)
        
        header_dict= mimage.get_attr_dict()
        header_dict["ptcl_source_coord_id"]= idx
        header_dict["data_path"]= sp_utilities.makerelpath(mpi_diff_bdb_obj.eman2db_dir, mpi_difference_stack)
        
        if stats_doc or do_norm:
            header_dict['old_mean'] = st[0]
            header_dict['old_sigma'] = st[1]
            
        new_bdb_dict[idx]= header_dict
        
        # Validate
        assert mimage.get_attr("ptcl_source_coord") == simage.get_attr("ptcl_source_coord"), "ERROR!! Mismatch %s != %s" % (
            mimage.get_attr("ptcl_source_coord"), simage.get_attr("ptcl_source_coord") )
        assert mimage.get_attr("ptcl_source_image") == simage.get_attr("ptcl_source_image"), "ERROR!! Mismatch %s != %s" % (
            mimage.get_attr("ptcl_source_image"), simage.get_attr("ptcl_source_image") )
    # End particle loop
    
    EMAN2db.db_close_dict(mpi_diff_bdb_obj.bdb_name)
    
    dimages= EMAN2.EMUtil.get_image_count(mpi_diff_bdb_obj.bdb_name)
    if dimages != nimages:
        mesg= "ERROR! %s should have %s images but instead has only %s. Exiting..." % (
                    mpi_diff_bdb_obj.bdb_name, nimages, dimages)
        print_log_msg(mesg, log)
        safe_exit()
    
    if RUNNING_UNDER_MPI and is_main and verbosity>=1: #### my_mpi_proc_id==main_mpi_proc 
        print_log_msg("Waiting for parallel processes to finish...\n", log)
    
    quick_barrier()
    if verbosity>=1 : print_log_msg("Finished subtracting reprojections\n", log)
    
    combined_diffimg_bdb_obj= BdbNames(combined_diffimg_bdb_stem, outdir)
    print('combined_diffimg_bdb_obj', combined_diffimg_bdb_obj.bdb_name)  #### DIAGNOSTIC
    
    # Write combined BDB (main process only)
    if RUNNING_UNDER_MPI and is_main: #### my_mpi_proc_id==main_mpi_proc:
        if verbosity>=1 : print_log_msg("Combining signal-subtracted stack files...", log)
        
        # Combine stacks
        num_merged_parts= merge_bdbs(
            outdir, 
            mpi_diff_bdb_stem, 
            num_mpi_procs, 
            combined_diffimg_bdb_obj, 
            verbosity=verbosity, 
            log=log
            )
            
        assert num_init_parts == num_merged_parts, "ERROR!! Mismatch input particles %s != %s merged particles" % (num_init_parts, num_merged_parts)
        
    # Save stats
    if is_main: #### my_mpi_proc_id == main_mpi_proc:
        if stats_doc or do_norm:
            if verbosity>=1 : print_log_msg("Reading mean values...", log)
            avg_list= EMAN2.EMUtil.get_all_attributes(combined_diffimg_bdb_obj.bdb_name, "old_mean")
            
            if verbosity>=1 : print_log_msg("Reading sigmas...", log)
            sigma_list= EMAN2.EMUtil.get_all_attributes(combined_diffimg_bdb_obj.bdb_name, "old_sigma")
            
            # Combine lists
            stats_tuples= list( zip(avg_list, sigma_list) )  # outer list will be required in Python 3
            
            # Convert to list of lists
            stats_list= [list(elem) for elem in stats_tuples]
            
            # Write to file
            sp_utilities.write_text_row(stats_list, stats_doc)
            
            if verbosity>=1 : print_log_msg("Wrote average and sigma before normalization to '%s'" % stats_doc, log)
        # End stats IF-THEN
    # End main IF-THEN
        
    quick_barrier()
    
def merge_bdbs(input_bdb_dir, input_stem_template, num_files, output_bdb_obj, verbosity=0, log=None):
    """
    Merges a series of BDBs. 
    Adapted from sp_tilt_import.py, which was adapted from sp_separate_class.py
    Assumes BDBs are consecutively numbered
    
    Modified 2021-01-05
    (TODO: Validate size, e.g., against another file)
    
    Arguments:
        input_bdb_dir : directory of input BDBs (not including EmAN2DB)
        input_stem_template : BDB stem (w/o leading 'bdb:') for input projection stacks
        num_files : number of BDBs to merge
        output_bdb_obj : instance of class BdbNames of combined stack
        verbosity :  How much information to write to screen (0..3)
        
    Returns:
        img_counter : number of entries merged
    """

    # Remove existing output BDB (will otherwise merge with pre-existing file)
    if os.path.exists(output_bdb_obj.bdb_path) : os.remove(output_bdb_obj.bdb_path)

    # We'll need this variable when updating the data_path in the header
    out_eman2db_dir= os.path.abspath(output_bdb_obj.eman2db_dir)
    
    # Open new database
    new_bdb_dict= EMAN2db.db_open_dict(output_bdb_obj.bdb_name)

    img_counter= 0  # initialize
    
    # Turn off progress bar if we're already printing information every iteration
    disableTF= verbosity!=2
    
    # Loop through BDBs
    for idx in tqdm.tqdm(range(num_files), unit='stack', disable=disableTF, file=sys.stdout):
        input_bdb_stem= input_stem_template.format(idx)
        input_bdb_obj= BdbNames(input_bdb_stem, input_bdb_dir)
        input_bdb_name= input_bdb_obj.bdb_name
        
        input_bdb_dict= EMAN2db.db_open_dict(input_bdb_name, ro=True)  
        
        # Get number of particles
        bdb_numparts= len(input_bdb_dict)
        if verbosity>=3 : print_log_msg("  Opened %s images from '%s'" % (bdb_numparts, input_bdb_name), log)
        
        # Loop through particles
        for partnum in range(bdb_numparts):
            # Read image header
            img_header= input_bdb_dict.get_header(partnum)
            
            # Update relative path between new BDB and image stack
            if 'data_path' in img_header:
                rel_path= update_data_path(
                        img_header['data_path'], 
                        os.path.join(input_bdb_dir, 'EMAN2DB'), 
                        out_eman2db_dir)
                img_header['data_path']= rel_path
            
            # (When are the following used?)
            img_header['data_source']= input_bdb_name
            img_header['data_n']= partnum
            
            # Write in new database
            new_bdb_dict[img_counter]= img_header
            img_counter += 1
        # End particle loop
        
        # Clean up
        EMAN2db.db_close_dict(input_bdb_name)
        del input_bdb_obj
    # End BDB loop
    
    # Close new database
    EMAN2db.db_close_dict(output_bdb_obj.bdb_name)
    print('1534: EMAN2db.db_close_dict', output_bdb_obj.bdb_name, type(output_bdb_obj.bdb_name) )  #### DIAGNOSTIC
    numoutimgs= EMAN2.EMUtil.get_image_count(output_bdb_obj.bdb_name)
    print('1536: Ignore warning "a bytes-like object is required, not NoneType", is an issue with EMAN2db.db_close_dict')  #### DIAGNOSTIC 
    assert numoutimgs == img_counter, "Uh oh!! merge_bdbs: %s != %s" % (numoutimgs, img_counter)
    if verbosity>=1 : 
        print_log_msg("Wrote %s images to '%s' from %s files\n" % (numoutimgs, output_bdb_obj.bdb_name, num_files), log)
        
    return img_counter
    
def update_data_path(curr_path, in_eman2db_dir, out_eman2db_dir):
    """
    Update relative data path between BDB and particle stack. 
    Copied from sp_separate_class.py
    
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
        print("ERROR!! update_data_path: %s does not exist!" % abs_img_path)
        printvars(['in_eman2db_dir', 'curr_path', 'out_eman2db_dir', 'abs_img_path'], quitTF=True)
    
    # Get relative path between source (BDB) and target (image stack)
    rel_path= sp_utilities.makerelpath(out_eman2db_dir, abs_img_path)
    
    return rel_path

def write_example_montages(orig_stack, proj_stack, diff_stack, numtestimgs, 
                           outdir, outstack, log=None, verbosity=0):
    """
    Writes montage of comparison images.
    
    Arguments:
        orig_stack : original BDB
        proj_stack : projection BDB
        diff_stack : difference BDB
        numtestimgs : number of test montages to write
        outdir : output directory
        outstack : output stack
        log : instance of Logger class
        verbosity : how much information to write to screen (0..3)
    """
    
    ##### Set up MPI
    ####if "OMPI_COMM_WORLD_SIZE" in os.environ:
        ####my_mpi_proc_id= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    ####else:
        ####my_mpi_proc_id= 0
        
    ####orig_stack= mpi_orig_template.format(my_mpi_proc_id)
    ####proj_stack= 'bdb:' + outdir + '#' + mpi_proj_bdb_stem.format(my_mpi_proc_id)
    ####diff_stack= 'bdb:' + outdir + '#' + mpi_diff_bdb_stem.format(my_mpi_proc_id)
    out_path= os.path.join(outdir, outstack)
    
    # Get number of images, and sanity check
    num_orig= EMAN2.EMUtil.get_image_count(orig_stack)
    num_proj= EMAN2.EMUtil.get_image_count(proj_stack)
    num_diff= EMAN2.EMUtil.get_image_count(diff_stack)
    assert num_orig==num_proj and num_proj==num_diff, "Uh oh! %s %s %s" % (num_orig,num_proj,num_diff)
    
    if numtestimgs>num_orig : numtestimgs= num_orig
    
    # Loop through images
    for idx in range(numtestimgs):
        # Read images
        orig_img= sp_utilities.get_im(orig_stack, idx)
        proj_img= sp_utilities.get_im(proj_stack, idx)
        diff_img= sp_utilities.get_im(diff_stack, idx)
        
        # Stitch montage
        montage_img= montage_scale([orig_img, proj_img, diff_img], scale=True)  # moved to sp_utilities in v1.3
        
        # Write to disk
        montage_img.write_image(out_path, idx)
    # End image loop
    
    if verbosity>=1:
        print_log_msg("Wrote %s example images to '%s'\n" % (numtestimgs, out_path), log)

def montage_scale(input_stack_or_list, ncol=None, marginwidth=0, bkgd=0, outfile=None, scale=False):
    """
    Generates montage of images into one image.
    Adapted from sxmontage.py and SPIDER's MN S
    Modified from in sp_utilities in SPHIRE 1.3
    2021-01-04 -- Updated for Python 3
    2021-01-05 -- Modified
    
    Arguments:
        input_stack_or_list : Stack of input images to merge into montage
        ncol : Number of images per row (default: all on one row)
        marginwidth : Margin width, pixels
        bkgd : Background value of montage
        outfile : Optional output file with montage output
        scale : Normalize images from 0..2
    Returns:
        montage : EMData object of image montage
    """
    
    if isinstance(input_stack_or_list, str): 
        input_list= EMAN2.EMData.read_images(input_stack_or_list)
    else:
        assert isinstance(input_stack_or_list, list), "Uh oh! Don't know type %s" % type(input_stack_or_list)
        input_list= input_stack_or_list
    
    # Get single-image dimensions
    nx= input_list[0].get_xsize()
    ny= input_list[0].get_ysize()
    
    # Get number of images and calculate montage dimensions
    numimgs= len(input_list)
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
            currentx= input_list[imgnum].get_xsize()
            assert currentx == nx, "ERROR!! Mismatched dimensions %s != %s" % (currentx, nx)
            [avg,var,fmin,fmax]= EMAN2.Util.infomask(input_list[imgnum], mask, True)
            
            img_norm= (input_list[imgnum] - fmin) * old_div(2., fmax - fmin)
        else:
            img_norm= input_list[imgnum]
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

def clean_projs(combined_proj_bdb_name, mpi_proj_bdb_stem, mpi_stack_template, outdir, mpioptions, 
                inmem=False, log=None, verbosity=0):
    """
    Deletes intermediate projection stacks.
    
    Arguments:
        combined_proj_bdb_name : Combined projection BDB
        mpi_proj_bdb_stem : BDB stem for each MPI process
        mpi_stack_template : stack template
        outdir : output directory
        mpioptions : (dict) parallelization options
        inmem : (boolean) compute projection in RAM
        log : instance of Logger class
        verbosity : how much information to write to screen (0..3)
    """
    
    # Loop through parallel processes
    for curr_proc_id in range(mpioptions['size']):
        # Get full path of '.bdb' file
        if mpioptions['use']:
            proj_bdb_obj= BdbNames(mpi_proj_bdb_stem.format(curr_proc_id), outdir)
        else:
            proj_bdb_obj= BdbNames(combined_proj_bdb_name)
        proj_bdb_path= proj_bdb_obj.bdb_path
        
        # If in memory, then there will be an MRC stack
        if inmem:
            # If not parallel, then no need to substitute process ID
            if mpioptions['use']:
                proj_stack_file= mpi_stack_template.format(curr_proc_id)
            else:
                proj_stack_file= mpi_stack_template
        # If not in memory, then we need to find the dimensions
        else:
            # If first process, get image dimensions (needed if not inmem)
            if curr_proc_id==0:
                img= sp_utilities.get_image(proj_bdb_obj.bdb_name, 0)
                xdim= img['nx']
                ydim= img['ny']
                zdim= img['nz']
            proj_stack_file= os.path.splitext(proj_bdb_path)[0] + "_%sx%sx%s" % (xdim,ydim,zdim)
        # End inmem if-then
        
        os.remove(proj_bdb_path)
        os.remove(proj_stack_file)
        if verbosity>=3:
            print_log_msg("  Deleted %s & %s" % (proj_bdb_path, proj_stack_file), log)
    # End process loop
    
def centershift(input1, input2, diffimgs, uncentered_file, inverted_file, volradius, init_centered_file, inparams, params_cp, params_cent, log, verbose, is_main):
    """
    Computes 3D reconstruction, centers, applies to image stack, and re-computes 3D reconstruction.
    
    Arguments:
        input1 : input map #1
        input2 : input map #2
        uncentered_file : output uncentered map
        options : (list) command-line options, primarily for sp_recons3d_n.py
        inverted_file : output contrast-inverted reconstruction, to make sure that there are no significant negative densities
        volradius : (float) radius of Gaussian blob to align reconstruction, pixels
        init_centered_file : output reconstruction after alignment to Gaussian blob
        params_cp : output initial alignment parameters
        params_cent : output alignment parameters after application of shifts required for centering
        log : instance of Logger class
        verbose : (boolean) -- whether to write to screen
        is_main : (boolean) if using multiple cores, some tasks only need to be performed once, not by all cores
    """
    
    uncentered_vol= average2(input1, input2, uncentered_file, log, verbose, is_main)
    
    # Negate input map
    volinvert= -1*uncentered_vol
    volinvert.write_image(inverted_file)
    if verbose: 
        mesg= "Wrote negated map '%s'. Check to make sure no density remains in subtracted region." % inverted_file
        print_log_msg(mesg, log)
    
    # Center map
    blob= make_blob(uncentered_vol, volradius, outfile=None)  # gets dimensions from input volume
    bestshift= align_vols(uncentered_vol, blob, numpeaks=2, log=log, verbose=verbose, is_main=is_main)
    uncentered_vol.transform(bestshift)
    uncentered_vol.write_image(init_centered_file)
    xshift= bestshift.get_params('spider').get('tx')
    yshift= bestshift.get_params('spider').get('ty')
    zshift= bestshift.get_params('spider').get('tz')
    if verbose: 
        print_log_msg("Shifted by (%s, %s, %s) to '%s'\n" % (xshift, yshift, zshift, init_centered_file), log)
    
    if inparams:
        # Copy alignment parameters
        if verbose : print_log_msg("Copying alignment parameters from '%s'" % inparams, log)
        shutil.copyfile(inparams, params_cp)
        
        # Combine shift with translation parameters
        if verbose: 
            print_log_msg("Applying shifts to alignment parameters '%s'" % params_cp, log)
            cmd=  "sp_process.py"
            args= "--transformparams=0,0,0,%s,%s,%s %s %s" % (xshift, yshift, zshift, params_cp, params_cent)
            print_log_msg(cmd + " " + args, log)
        
        # Adapted from sp_process.py
        python_version, devnull= setup_subprocess()
        system_call_23(cmd, args, stdout=devnull)
        transf= [0.0, 0.0, 0.0, float(xshift), float(yshift), float(zshift)]
        sp_utilities.write_text_row(sp_utilities.rotate_shift_params(sp_utilities.read_text_row(params_cp), transf), params_cent)
        if verbose : print_log_msg("Finished applying shifts, recorded in '%s'\n" % params_cent, log)
        
    else:  # if no input parameters specified
        if diffimgs:
            # Export current parameters for data stack
            if verbose:
                print_log_msg("Exporting alignment parameters %s from stack %s" % (params_cp, diffimgs), log)
                cmd= "sp_header.py %s --params=xform.projection --export=%s" % (diffimgs, params_cp)
                print_log_msg(cmd, log)
            sp_applications.header(diffimgs, 'xform.projection', fexport=params_cp)
        if verbose : print_log_msg("Wrote alignment parameters to stack '%s'\n" % diffimgs, log)
            
    if diffimgs:
        # Read in centered parameters into headers of the new stack
        if verbose:
            print_log_msg("Importing alignment parameters '%s'" % params_cent, log)
            cmd= "sp_header.py %s --params=xform.projection --import=%s" % (diffimgs, params_cent)
            print_log_msg(cmd, log)
        sp_applications.header(diffimgs, 'xform.projection', fimport=params_cent)
        if verbose : print_log_msg("Wrote alignment parameters to stack '%s'\n" % diffimgs, log)
    
def make_blob(dimensions, radiuspx, outfile=None):
    """
    Generates a Gaussian blob.
    
    Arguments:
        dimensions : (list) list of 1 or 3 elements of volume dimension, or EMData object from which dimensions will be copied
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
    
    cmd= "Finding shift parameters...\n" + \
        " "*TIMESTAMP_LENGTH + "   PEAK  XPEAK  YPEAK  ZPEAK  RPEAK  XORIG  YORIG  ZORIG"
    print_log_msg(cmd, log)
    
    for line in peaks: print_log_msg("".join(format(item, "7.1f") for item in line), log)

    xshift= -peaks[0][5]
    yshift= -peaks[0][6]
    zshift= -peaks[0][7]
    shift= EMAN2.Transform({"type":"spider","phi":0,"theta":0,"psi":0,"tx":xshift,"ty":yshift,"tz":zshift})
    
    return shift

def parse_command_line():
    """
    Parse the command line.  Adapted from sp_mask.py

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
        'mode', 
        type=str, 
        help='mode (options: avgfilt, sp_mask, projsubtract, centershift)')
    
    parser.add_argument(
        '--verbosity', "-v", 
        type=int, 
        default=2, 
        help='Increase verbosity (0..3)')
    """
    Verbosity levels:
        0: Only warnings
        1: Basic
        2: Progress bar over classes
        3: Merged classes
    """
    
    parser.add_argument(
        '--outdir', 
        type=str, default='.', 
        help='Output directory (default: .)')
    
    parser.add_argument(
        '--apix', type=float, 
        default=None, 
        help='Pixel size, Angstroms')
    
    group_avgfilt= parser.add_argument_group(
        title='AVGFILT options (mode avgfilt)',
        description='Will average two maps (if both available) and optionally low-pass filter result.')
    
    group_avgfilt.add_argument(
        '--avol1', 
        type=str, 
        default=None, 
        help='Input volume #1')
    
    group_avgfilt.add_argument(
        '--avol2', 
        type=str, 
        default=None, 
        help='Input volume #2')
    
    group_avgfilt.add_argument(
        '--filtrad', 
        type=float, 
        default=None, 
        help='Radius for low-pass filtration, px^-1 or Angstroms (if apix specified)')
    
    group_softmask= parser.add_argument_group(
        title='Adaptive-masking options (mode sp_mask)',
        description='Mask will be generated by alternately low-pass filtering and thresholding input map.')
    
    group_softmask.add_argument(
        '--map2mask', 
        type=str, 
        default=None, 
        help='Map out of which to make mask')
    
    group_softmask.add_argument(
        '--fullmap', 
        type=str, 
        default=None, 
        help='Map to which to apply mask')
    
    group_softmask.add_argument(
        '--mapthresh', 
        type=float, 
        default=None, 
        help='Map threshold for masking (default: +1 s.d.)')
    
    group_softmask.add_argument(
        '--maskthresh', 
        type=float, 
        default=None, 
        help='Map threshold for masking (0..1)')
    
    group_softmask.add_argument(
        '--ndilation', 
        type=int, 
        default=3, 
        help='Number of dilations')
    
    group_softmask.add_argument(
        '--edge_width', 
        type=int, 
        default=5, 
        help='Number of dilations')
    
    group_softmask.add_argument(
        '--masklowpass', 
        type=int, 
        default=5, 
        help='Low-pass filter radius, Angstroms (or if, apix not provided, pixels^-1)')
    
    group_softmask.add_argument(
        '--falloff', 
        type=float, 
        default=0.01, 
        help='Low-pass filter falloff in absolute frequency')
    
    group_softmask.add_argument(
        '--allow_disconnected', '--ad',
        action='store_true',
        default=False,
        help='Allow disconnected regions in the mask.')
    
    group_projsubtract= parser.add_argument_group(
        title='Projection-subtraction options (mode projsubtract)',
        description='To subtract projections of masked map from experimental images.')
    
    group_projsubtract.add_argument(
        '--origparts', 
        type=str, 
        default=None, 
        help='Original particle stack')
    
    group_projsubtract.add_argument(
        '--map2subtract', 
        type=str, 
        default=None, 
        help='Map to subtract')
    
    group_projsubtract.add_argument(
        '--projparams', 
        type=str, 
        default=None, 
        help='Projection parameter text file')

    group_projsubtract.add_argument(
        '--nonorm', 
        action="store_true", 
        help='Turn off normalization of images')

    group_projsubtract.add_argument(
        '--nmontage', 
        type=int, 
        default=0, 
        help='Write this number of test images of projections and subtracted images')

    group_projsubtract.add_argument(
        '--inmem', 
        action="store_true", 
        help='Store reprojections in memory (may be large)')

    group_projsubtract.add_argument(
        '--saveprojs', 
        action="store_true", 
        help='Save reprojections')

    group_projsubtract.add_argument(
        '--stats', 
        action="store_true", 
        help='Write average and sigma for input particles')

    group_projsubtract.add_argument(
        '--normalize', 
        action="store_true", 
        help='DEPRECATED: Normalize subtracted images (is now on by default)')

    group_centershift= parser.add_argument_group(
        title='CENTERSHIFT options (mode centershift)',
        description='To center reconstruction.')
    
    group_centershift.add_argument(
        '--cvol1', 
        type=str, 
        default=None, 
        help='Input volume #1')
    
    group_centershift.add_argument(
        '--cvol2', 
        type=str, 
        default=None, 
        help='Input volume #2')
    
    group_centershift.add_argument(
        '--diffimgs', 
        type=str, 
        default=None, 
        help='Difference image stack')
    
    group_centershift.add_argument(
        '--shiftparams', 
        type=str, 
        default=None, 
        help='Refinement parameter text file')
    
    group_centershift.add_argument(
        '--volradius', 
        type=float, 
        default=None, 
        help='Radius of map, pixels or (if apix specified) Angstroms')
    
    group_centershift.add_argument(
        "--prefix", 
        type=str, 
        default="vol",  
        help="Prefix for output reconstructions (default: vol)" )
    
    return parser.parse_args()

def main():
    options= parse_command_line()

    #print('options', options)
    #exit()
    
    RUNNING_UNDER_MPI= "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
            mpi.mpi_init( 0, [] )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp( "Start" )
    
    signalsubtract(options, outdir=options.outdir, verbosity=options.verbosity)

    sp_global_def.print_timestamp( "Finish" )
    if RUNNING_UNDER_MPI:
            mpi.mpi_finalize()

if __name__ == "__main__":
    main()
