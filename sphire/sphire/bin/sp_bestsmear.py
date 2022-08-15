#!/usr/bin/env python
from __future__ import print_function

import os
import mpi
import argparse
import inspect
import six
import subprocess
import sys
import numpy 
import re 
import json
import datetime
import random
import inspect
import six
import EMAN2
from sphire.libpy import sp_global_def 
from sphire.libpy import sp_logger
from sphire.libpy import sp_utilities
from matplotlib import pyplot as plt

USAGE= """ 
PURPOSE:
Keeps the particles with the lowest smear number (highest certainty).

General usage:
%s <meridien_parameters> <output_directory>
General, required command-line parameters:
  Meridien alignment parameters, e.g., final_params_027.txt
Outputs:
  smear_info.txt : Smear number for each particle
  plot_smear.png : Plot of lowest of highest smear numbers

General options:
%s <meridien_parameters> <output_directory> --iter=<iteration_num> --threads=<num_threads> --verbosity=<verbosity>
Parameters:
  --iter : Iteration number for Meridien parameter file (default: automatic)
  --threads : Number of threads used in Meridien (default: automatic)
  --verbosity : Increase verbosity (range: 0-5)

To apply a threshold smear number:
%s <meridien_parameters> <output_directory> --max_smear=<maximum_smear>
Parameter:
  --max_smear : Maximum smear to permit

To apply a threshold smear number:
%s <meridien_parameters> <output_directory> --max_smear=<maximum_smear> --stack <input_stack>
Parameter:
  --stack : Stack from which subset with the lowest smear will be retained

""" % ((__file__,)*4)

MODIFIED= "Modified 2021-02-19"
TIMESTAMP_LENGTH= 23

# Meridien inputs
CHUNK_FILE= 'main%03d/chunk_%01d_%03d.txt'
JSON_PARAMS= 'main%03d/oldparamstructure/oldparamstructure_%01d_%03d_%03d.json'

# Outputs
SMEAR_DOC= 'smear_info.txt'
PLOT_FILE= 'plot_smear.png'
GOOD_PARTS= 'good_parts.txt'
GOOD_STACK= 'bestsmear'

def setup_mpi(usempi):
    """
    Set up MPI.
    
    Argument:
        usempi : (boolean) Whether to use Message Passing Interface
        
    Returns:
        mpidict : dictionary with keys:
            'use' : (boolean) whether to use MPI
            'myid' : MPI process number
            'main_node' : main MPI process number
            'is_main' : (boolean) whether current process is main
    """
    
    if usempi:
        number_of_proc= mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        myid= mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        
        #  choose a random node as main
        main_node= 0
        if(myid == 0) : main_node= random.randint(0,number_of_proc-1)
        main_node= mpi.mpi_bcast(main_node, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)	
    else:
        myid= main_node= 0
        # (IF-THEN statements will be simpler if I pretend non-MPI runs are the main node.)
    
    is_main= myid == int(main_node)  # boolean: whether current process is the main one
    
    mpidict= {'use':usempi, 'myid':myid, 'main_node':main_node, 'is_main':is_main}
    return mpidict
    
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
        print("WARNING: Using old logger.py library")
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

def print_log_msg(msg, log=None, verbose=False, main=True):
    """
    Prints messages to log file and, optionally, to the screen.
    
    Arguments:
        msg : Message to write
        log : Instance of Logger class
        verbose : (boolean) Whether to write to screen
        is_main : (boolean) If using MPI, some tasks only need to be performed once, not by all cores
    """
    
    if main:
        if verbose: print(msg)
        if log: 
            log.add(msg)
        else:
            print(msg)

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
        msg= "%s: '%s'" % (variable, dictionary[variable])
        if typeTF: msg+= " %s" % type(dictionary[variable])
        print(msg)
        
    del frame
    
    if quitTF: 
        print('\nExiting printvars...')  # reminder in case I forget to take out the quit flag
        exit()

def system_call_23(cmd, args, lenient=False, stdout=None, stderr=None, log=None, verbose=False):
    """
    Runs subprocess safely.
    
    Arguments:
        cmd : Executable
        args : Command-line arguments
        lenient : (boolean) Will simply print warning in case of error
        stdout : Where to direct standard out
        stderr : Where to direct standard error
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
    
    try:
        if version_float < 3.5:
            printvars('version_float')
            subprocess.check_call([cmd, args], stdout=stdout, stderr=stderr, shell=True)
            # (shell=True is risky, but I couldn't get it to work without)
        else:
            printvars('version_float')
            subprocess.run([cmd] + args.split(), stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError:
        if not lenient:
            print("ERROR!! Cannot execute '%s'\n" % cmdline)
            exit()
        else:
            print("WARNING! Cannot execute '%s'\n" % cmdline)

def best_smear(params_doc, outdir='.', max_smear=None, verbosity=0, options=None, usempi=False):
    """
    FUNCTION DESCRIPTION.
    
    Arguments:
        outdir : Output directory
        options : (list) Command-line options
        verbosity : How much information to write to screen
        usempi : (boolean) Whether to use Message Passing Interface
    """
    
    #Set up MPI
    mpioptions= setup_mpi(usempi)
    
    # Set output directory and log file name
    prepare_outdir(outdir, main=mpioptions['is_main'])
    sp_global_def.write_command(outdir)
    log, _= prepare_log(outdir)
    
    # Wait
    if usempi: mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    
    """
    Read the number of smeared references for each image.
    Adapted from sp_bestim.py

    In the Meridien/main###/oldparamstructure directory, 
        there is a file for each MPI thread called oldparamstructure_0_???_###.json, 
        where ??? is the thread number.
    In the JSON file, each particle is assigned a list which contains:
        [ particle#, [0.0], [[data_reference1], [data_reference2], ...] ]
    The number of projections to be smeared is what we are looking for in the first pass.
    """

    # Define filenames
    refine_dir= os.path.dirname(params_doc)
    chunk_template= os.path.join(refine_dir, CHUNK_FILE)
    oldparam_template= os.path.join(refine_dir, JSON_PARAMS)
    smeardoc= os.path.join(outdir, SMEAR_DOC)
    plotfile= os.path.join(outdir, PLOT_FILE)
    goodsel= os.path.join(outdir, GOOD_PARTS)
    newstack= 'bdb:' + outdir + "#" + GOOD_STACK

    # Get iteration number
    if not options.iter:
        iter_num= get_filenumber(options.inparams)
    else:
        iter_num= options.iter
    
    # Get number of iterations
    if not options.threads:
        num_threads= get_numthreads(oldparam_template, iter_num)
    else:
        num_threads= options.threads
    
    if verbosity >=1 : 
        print_log_msg("Verbosity level set to %s (max 5)" % verbosity, log)
        print_log_msg('Reading refinement data from directory: %s' % refine_dir, log)
        
        if max_smear:
            print_log_msg('Excluding particles with a smear number greater than %s' % max_smear, log)
        
        if options.stack:
            if max_smear: 
                print_log_msg("Will write subset of '%s'" % options.stack, log)
            else:
                print_log_msg("WARNING! Option '--stack' does nothing without '--max_smear' ", log)
        
        if options.iter:
            print_log_msg("Using Meridien iteration #%s" % iter_num, log)
        else:
            print_log_msg("Automatically using Meridien iteration #%s" % iter_num, log)
        
        if options.threads:
            print_log_msg("Assuming %s parallel threads" % num_threads, log)
        else:
            print_log_msg("Found %s parallel threads" % num_threads, log)
        
        print()
        print_log_msg('Getting smear information...', log)
    
    # Get number of total particles
    num_globparts= len( numpy.genfromtxt(params_doc) )

    # Initialize list with smear number for each particle
    smear_list= [None]*num_globparts

    # Initialize global particle counter
    globpart_counter= 0

    # Loop through chunks
    for chunk_num in range(2):
        # Read list of particles for this chunk
        chunk_partfile= chunk_template % (iter_num, chunk_num, iter_num)
        chunk_partlist= numpy.genfromtxt(chunk_partfile)
        if verbosity >= 3 : print_log_msg('  %s: %s parts' % (chunk_partfile, len(chunk_partlist) ), log)

        # Initialize particle-counter
        chunk_partnum= -1

        # Loop through threads
        for thread_num in range(num_threads):
            # Open JSON parameter file
            oldparam_file= oldparam_template % (iter_num, chunk_num, thread_num, iter_num)
            with open(oldparam_file) as read:
                num_thread_parts= json.load(read)
            
            # Loop through particle data
            for part_data in num_thread_parts:
                thread_partnum= part_data[0]
                
                # Get the number of references included in the smearing
                num_smear= len(part_data[2])
                # Increment global particle counter
                chunk_partnum += 1
                
                # Read global particle number
                glob_partnum= int(chunk_partlist[chunk_partnum])
                
                if verbosity >= 5 : 
                    msg= '  thread# %s, thread_part# %s, chunk_part# %s, glob_part# %s, #smear %s' % (
                        thread_num, thread_partnum, chunk_partnum, glob_partnum, num_smear 
                        )
                    print_log_msg(msg, log)
                
                smear_list[glob_partnum]= num_smear
            # End particle-loop
            
            if verbosity >= 4: 
                print_log_msg('  %s, %s particles' % (oldparam_file, len(num_thread_parts) ), log)
        # End thread-loop

        chunk_partnum += 1
        globpart_counter += chunk_partnum
        if verbosity >= 3 :
            mesg= '  Chunk #%s: Got smear information for %s particles\n' % (chunk_num, chunk_partnum)
            print_log_msg(mesg, log)
    # End chunk-loop

    #write_text_row(zip(range(num_globparts), smear_list), smeardoc)  
    sp_utilities.write_text_row(smear_list, smeardoc) 
    # (TODO: I'd like to include the particle# also, but didn't know how. -Tapu)
    
    if verbosity >= 1: 
        mesg= "Wrote smear information for %s particles to '%s'" % (globpart_counter, smeardoc)
        print_log_msg(mesg, log)
    
    if num_globparts != globpart_counter: 
        sp_global_def.ERROR(
            "\nFound %s particles in %s, but %s particles in %s\n" % 
            (num_globparts, params_doc, globpart_counter, os.path.dirname(oldparam_file) ), 
            __file__, 1)
        exit()
    
    sort_smear= sorted(smear_list)
    smear_ranks= [x+1 for x in range( len(smear_list) )]
    
    try:
        if max_smear != None : plt.subplot(2,1,1)
        plt.scatter(smear_ranks, sort_smear, s=1)
        if not max_smear : plt.xlabel('Particle rank')
        plt.ylabel('Smear number')
        if max_smear : 
            plt.title('Smear number before & after cutoff')
        else:
            plt.title('Smear number, no cutoff')
            plt.savefig(plotfile)
            if verbosity >=1:
                print_log_msg('Wrote plot of smear number without cutoff to %s\n' % plotfile, log)

    except RuntimeError:
        print_log_msg("WARNING! Can't create plot (a machine-specific problem). Skipping...", log)

    # Remove high-smear particles
    if max_smear:
        # Initialize lists
        keep_list= []
        kept_smear= []
        
        # Save only particles not exceeding maximum specified smear
        for idx in range( len(smear_list) ):
            if smear_list[idx] <= max_smear:
                keep_list.append(idx)
                kept_smear.append(smear_list[idx])
                
        sp_utilities.write_text_row(keep_list, goodsel)
        
        if verbosity >=1:
            mesg= "Wrote %s out of %s particles not exceeding smear of %s to '%s'" % (
                len(keep_list), len(smear_list), max_smear, goodsel
                )
            print_log_msg(mesg, log)
        
        sort_smear= sorted(kept_smear)
        smear_ranks= [x+1 for x in range( len(keep_list) )]
    
        try:
            plt.subplot(2,1,2)
            plt.scatter(smear_ranks, sort_smear, s=1)
            plt.xlabel('Particle rank')
            plt.ylabel('Smear number')
            plt.savefig(plotfile)
            if verbosity>=1:
                print_log_msg("Wrote plot of smear number before & after cutoff to '%s'\n" % plotfile, log)
        except RuntimeError:
            pass  # already wrote warning above
    
        if options.stack:
            # Write substack
            cmd= "e2bdb.py %s --makevstack %s --list %s" % (options.stack, newstack, goodsel)
            if verbosity>=1: 
                print_log_msg("Writing substack of '%s' " % options.stack, log)
                print_log_msg(cmd, log)
            os.system(cmd)  # (TODO: use safer system_call_23)
            
            num_new_parts= EMAN2.EMUtil.get_image_count(newstack)
            if verbosity>=1: 
                print_log_msg("Wrote %s images to '%s'" % (num_new_parts, newstack), log)
    
    if verbosity>=1 : print_log_msg("Done!\n", log)
    ####if mpioptions['is_main']: print("Done!")
    
def get_filenumber(filename):
    """
    Extracts filenumber from filename.
    
    Argument:
        filename : Assumes single number per filename
    """
    
    # Separate filename from full path (which might have numbers)
    basename = os.path.basename(filename)
    
    # Find numbers
    filenumbers = re.findall(r'\d+', basename)
    # TODO: Check for nonstandard name/multiple numbers
    
    filenumber = int(filenumbers[0])
    
    return filenumber

def get_numthreads(template, iter_num):
    """
    Get number of threads by looking for files numbered by thread.

    Arguments:
            template : Filename template, of the form main%03d/oldparamstructure/oldparamstructure_%01d_%03d_%03d.json'
            iter_num : Iteration number in template
    """

    # Initialize counter
    thread_num = 0

    while True:
            # Substitute values into filename template
            file = template % (iter_num, 0, thread_num, iter_num)

            # Check if file exists
            check = os.path.exists(file)

            if not check: break

            thread_num += 1

    max_thread = thread_num

    if max_thread == 0:
            sp_global_def.ERROR("\nCouldn't find Meridien file %s\n" % file, __file__, 1)
            exit()

    return max_thread

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
        'inparams', 
        type=str, 
        help='Input Meridien parameter file')
    
    parser.add_argument(
        'outdir', 
        type=str, 
        default='.', 
        help='Output directory (default: .)')
    
    parser.add_argument(
        '--max_smear', '-m',
        type=int, 
        default=None,      
        help='Maximum smear number to allow (default: keep all)')
    
    parser.add_argument(
        '--stack',
        type=str, 
        default=None,      
        help='Stack for which to take subset')
    
    parser.add_argument(
        '--iter', '-i',
        type=int, 
        default=None,      
        help='Iteration number (default: extract from params filename)')
    
    parser.add_argument(
        '--threads', '-t',
        type=int, 
        default=None,      
        help='Number of threads in Meridien (default: extract automatically)')
    
    parser.add_argument(
        '--verbosity', "-v", 
        type=int, 
        default=3,      
        help='Increase verbosity (range 0-5)')
    
    return parser.parse_args()

def main():
    options= parse_command_line()
    
    ##print args, options  # (Everything is in options.)
    #print(options)
    ##print('LOGFILE',sp_global_def.LOGFILE)
    #exit()

    RUNNING_UNDER_MPI= "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
            mpi.mpi_init( 0, [] )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp( "Start" )	

    best_smear(
        options.inparams, 
        outdir=options.outdir, 
        max_smear=options.max_smear, 
        verbosity=options.verbosity,
        options=options
        )

    sp_global_def.print_timestamp( "Finish" )
    if RUNNING_UNDER_MPI:
            mpi.mpi_finalize()

if __name__ == "__main__":
    main()
