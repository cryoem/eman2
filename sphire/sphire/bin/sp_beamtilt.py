#!/usr/bin/env python
from __future__ import print_function

import os
import datetime
import mpi
import random
import argparse
import inspect
import six
import numpy as np
import EMAN2
import math
import re
import tqdm
import collections
import EMAN2db
import pathlib
import glob
import sys
import itertools
import subprocess
import scipy
import pandas

from sphire.libpy import sp_global_def
from sphire.libpy import sp_logger
from sphire.libpy import sp_applications
from sphire.libpy import sp_utilities

TIMESTAMP_LENGTH = 23
MAX_VERBOSE = 7
EXTRACT_DIR = 'FromSTAR'
UNTILT_DIR = 'Corrected'
COMBINEDBDB = 'all_corrected'

USAGE = """ 
Applies beam-tilt correction to particles.

  %s <output_directory> <options>
General, required command-line parameters:
  Output directory
Output:
  bdb:all_corrected : tilt-corrected image stack

General options:
  %s <output_directory> --input_stack=INPUT_STACK --star_file=STAR_FILE --verbosity --debug
Parameters:
  --star_file : Input STAR file
  --input_stack : Input particle stack
  --verbosity : Increase verbosity (0..7)
  --debug : Flag for dry run

To get beam-tilt parameters from a STAR file and apply correction to a BDB stack:
  %s <output_directory> --input_stack=INPUT_STACK --star_file=STAR_FILE

To provide beam-tilt parameters directly, without a STAR file:
  %s <output_directory> --input_stack=INPUT_STACK --kv=KV --cs=CS --apix=APIX --tiltx=TILTX --tilty=TILTY
Parameters:
  --kv : Accelerating voltage, kV
  --cs : Spherical aberration constant, mm
  --apix : Pixel size, Angstroms
  --tiltx : Beam tilt along x, milliradians
  --tilty : Beam tilt along y, milliradians

If particles need to be extracted from micrographs:
  %s <output_directory> --star_file=STAR_FILE --micdir=MICDIR  --box_size=BOX_SIZE--coordsdir=COORDSDIR --partdir=PARTDIR
Parameters:
  --micdir : Input directory where micrographs are located 
  --box_size : Dimension for extracted particles 
  --coordsdir : Output directory where particle coordinates will be written
  --partdir : Output directory where particles will be written

""" % ((__file__,) * 5)

MODIFIED = "Modified 2021-02-18"


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
        self.bdb_stem = bdb_stem
        if dir_above_eman2db != None:
            self.bdb_dir = dir_above_eman2db
        else:
            self.bdb_dir = '.'
        self.eman2db_dir = os.path.join(self.bdb_dir, 'EMAN2DB')
        self.bdb_name = self.stem2name()
        self.bdb_path = self.stem2path()

    def stem2name(self):
        if self.bdb_dir != None:
            name = 'bdb:' + self.bdb_dir + "#" + self.bdb_stem
        else:
            name = 'bdb:' + self.bdb_stem

        return name

    def stem2path(self):
        if self.bdb_dir != None:
            path = os.path.join(self.bdb_dir, 'EMAN2DB', self.bdb_stem + '.bdb')
        else:
            path = os.path.join('EMAN2DB', self.bdb_stem + '.bdb')

        return path

    def list_info(self):
        info = []
        info.append(self.bdb_stem)
        info.append(self.bdb_dir)
        info.append(self.eman2db_dir)
        info.append(self.bdb_name)
        info.append(self.bdb_path)

        return info


class StarFile:
    """
    Author: Markus Stabrin
    Adapted and annotated, 2020-07-01, TRS

    This script takes a relion 3.1 star file and converts it into a relion 3.0 compatible star file.

    Functions:
        __init__() : Initializes
        analyse_star_file(() : Determines line numbers for each block of data
        read_tag() : Populates 'imported_content' with data
        read_with_loop() :  Reads data when block starts with 'loop_'
        read_without_loop() : Reads data when block doesn't start with 'loop_'
        __getitem__() : Returns data
        __setitem__() : Imports data
        write_star() : Writes data to disk for specified tag(s)

    Variables:
        star_file : Nmae of STAR file
        imported_content : Data, as a pandas DataFrame
        line_dict : (dict) Line numbers for each section of each block

    """

    def __init__(self, star_file):
        self.star_file = star_file
        self.imported_content = {}
        self.line_dict = {}
        self.analyse_star_file()  # parses STAR file

    def analyse_star_file(self):
        """
        Populates self.line_dict with line numbers for each section of each block.

        line_dict : Dictionary whose keys are a block of the STAR file, with 'data_' removed (e.g., 'data_optics' -> 'optics')
            keys in each block:
                block : Line numbers for entire block
                header : Line numbers for header
                content : Line numbers for data
                is_loop : (boolean) Whether block starts with 'loop_'
        """

        with open(self.star_file) as read:
            content = read.read()

        # https://regex101.com/r/D7O06N/1
        for tag_match in re.finditer('^data_([^\s]*)\s*$', content, re.M):
            tag = tag_match.group(1)
            self.line_dict[tag] = {
                'block': [None, None],
                'header': [None, None],
                'content': [None, None],
                'is_loop': None,
            }

            current_flag = 0
            prev_content = content[:tag_match.start() + current_flag]
            current_content = content[tag_match.start() + current_flag:]
            current_flag += tag_match.start()

            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['block'][0] = \
                len(re.findall('\n', prev_content)) + 1

            # https://regex101.com/r/h7Wm8y/2
            header_match = re.search(
                '((?:(?:loop_\s*)?^_.*$\r?\n?)+)',
                current_content,
                re.M
            )
            prev_content = content[:header_match.start() + current_flag]
            current_content = content[header_match.start() + current_flag:]
            current_flag += header_match.start()

            self.line_dict[tag]['is_loop'] = header_match.group(1).startswith('loop_')
            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['header'][0] = \
                len(re.findall('\n', prev_content)) + 1 + self.line_dict[tag]['is_loop']

            prev_content = content[:header_match.end() + current_flag - header_match.start()]
            current_content = content[header_match.end() + current_flag - header_match.start():]
            current_flag += header_match.end() - header_match.start()
            # https://regex101.com/r/4o3dNy/1/
            self.line_dict[tag]['header'][1] = \
                len(re.findall('\n', prev_content))

            if not self.line_dict[tag]['is_loop']:
                self.line_dict[tag]['content'] = self.line_dict[tag]['header']
            else:
                self.line_dict[tag]['content'][0] = self.line_dict[tag]['header'][1] + 1
                # https://regex101.com/r/HYnKMl/1
                newline_match = re.search('^\s*$', current_content, re.M)

                prev_content = content[:newline_match.start() + current_flag]
                current_content = content[newline_match.start() + current_flag:]
                current_flag += newline_match.start()

                # https://regex101.com/r/4o3dNy/1/
                self.line_dict[tag]['content'][1] = \
                    len(re.findall('\n', prev_content))

            self.line_dict[tag]['block'][1] = self.line_dict[tag]['content'][1]

            self.read_tag(tag, self.line_dict[tag])

    def read_tag(self, tag, line_dict):
        """
        Populates self.imported_content with data.
        """

        if not line_dict['is_loop']:
            data = self.read_without_loop(line_dict)
        else:
            data = self.read_with_loop(line_dict)
        self.imported_content[tag] = data

    def read_with_loop(self, line_dict):
        """
        Reads data when block starts with 'loop_'.
        """

        header_names = pandas.read_csv(
            self.star_file,
            usecols=[0],
            skiprows=line_dict['header'][0] - 1,
            nrows=line_dict['header'][1] - line_dict['header'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
            squeeze=True,
        )

        return pandas.read_csv(
            self.star_file,
            index_col=None,
            names=header_names,
            skiprows=line_dict['content'][0] - 1,
            nrows=line_dict['content'][1] - line_dict['content'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
        )

    def read_without_loop(self, line_dict):
        """
        Reads data when block doesn't start with 'loop_'.
        """

        return pandas.read_csv(
            self.star_file,
            index_col=0,
            names=['', '0'],
            skiprows=line_dict['content'][0] - 1,
            nrows=line_dict['content'][1] - line_dict['content'][0] + 1,
            skip_blank_lines=False,
            header=None,
            delim_whitespace=True,
        ).transpose()

    def __getitem__(self, tag):
        return self.imported_content[tag]

    def __setitem__(self, data, tag):
        self.imported_content[tag] = data

    def write_star(self, star_file, tags):
        """
        Writes data to disk for specified tag(s).
        """

        for idx, tag in enumerate(tags):
            if idx == 0:
                mode = 'w'
            else:
                mode = 'a'
            df = self.imported_content[tag]
            is_loop = self.line_dict[tag]['is_loop']

            if is_loop:
                export_header = '\ndata_\n\nloop_\n' + '\n'.join([
                    '{} #{}'.format(entry, idx)
                    for idx, entry
                    in enumerate(df, 1)
                ])

                with open(star_file, mode) as write:
                    write.write(f'{export_header}\n')
                df.to_csv(star_file, sep='\t', header=False, index=False, mode='a')


def setup_mpi(usempi):
    """
    Set up MPI.
    Modified 2020-08-15

    Argument:
        usempi : (boolean) Whether to use Message Passing Interface

    Returns:
        mpidict : dictionary with keys:
            'use' : (boolean) whether to use MPI
            'myid' : MPI process number
            'main_node' : main MPI process number
            'is_main' : (boolean) whether current process is main
            'size' : number of processes
    """

    if usempi:

        mpi_size = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)

        #  choose a random node as main
        main_node = 0
        if (myid == 0): main_node = random.randint(0, mpi_size - 1)
        main_node = mpi.mpi_bcast(main_node, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)
    else:
        myid = main_node = 0
        # (IF-THEN statements will be simpler if I pretend non-MPI runs are the main node.)

        mpi_size = 1

    is_main = myid == int(main_node)  # boolean: whether current process is the main one

    mpidict = {'use': usempi, 'myid': myid, 'main_node': main_node, 'is_main': is_main, 'size': mpi_size}

    return mpidict


def quick_barrier():
    """
    Synchronizes parallel processes before continuing
    Safe for non-MPI calls also
    """

    if "OMPI_COMM_WORLD_SIZE" in os.environ: mpi.mpi_barrier(mpi.MPI_COMM_WORLD)


def safe_exit(verbosity=0):
    """
    Properly closes MPI before exiting, preventing ugly warnings
    Safe for non-MPI calls also

    Argument:
            verbosity : controls how much information to write to screen (0..3)
    """

    if "OMPI_COMM_WORLD_SIZE" in os.environ:
        my_rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        if verbosity >= 3:
            print("Process %s: Synchronizing with barrier" % my_rank)
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        if verbosity >= 2:
            print("Process %s: Barrier reached" % my_rank)
        mpi.mpi_finalize()
    else:
        my_rank = 0

    if verbosity >= 1:
        print("Process %s: Finalized" % my_rank)

    sp_global_def.print_timestamp("Finish")

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

    logname = "log_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + ".txt"
    logname = os.path.join(outdir, logname)

    # May be using old version of logger.py
    try:
        if verbose:
            log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), base_logger2=sp_logger.BaseLogger_Print(),
                                   file_name=logname)
            verbose = False  # logger output will be echoed to screen
        else:
            log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        print("WARNING: Using old logger.py library")
        log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())  # , file_name=logname)
        logname = 'log.txt'

    if main:
        print("Writing log file to %s" % logname)

        progbase = os.path.basename(__file__).split('.')[0].upper()
        length = len(progbase) + 4

        log.add("\n" +
                " " * TIMESTAMP_LENGTH + "*" * length + "\n" +
                " " * TIMESTAMP_LENGTH + "* " + progbase + " *\n" +
                " " * TIMESTAMP_LENGTH + "*" * length)

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

    import inspect
    import six

    if type(variables) is list:
        # Weird things happen if
        assert isinstance(variables[0], six.string_types), "UH OH!! Passed non-string %s instead of variable name" % \
                                                           variables[0]

        variable_list = variables
    elif isinstance(variables, six.string_types):  # six works with both Pythons 2 & 3
        variable_list = [variables]
    else:
        print("ERROR!! Don't know how to deal with type %s" % type(variables))
        exit()

    frame = inspect.currentframe()
    dictionary = frame.f_back.f_locals

    print("")
    for variable in variable_list:
        msg = "%s: '%s'" % (variable, dictionary[variable])
        if typeTF: msg += " %s" % type(dictionary[variable])
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
    python_version = sys.version.split()[0]

    # Keep up to one decimal
    version_float = float('.'.join(python_version.split('.', 2)[:2]))

    # Build command line
    cmdline = "%s %s" % (cmd, args)

    if verbose: print_log_msg(cmdline, log)

    try:
        if version_float < 3.5:
            subprocess.check_call([cmd, args], stdout=stdout, stderr=stderr, shell=True)
            # (shell=True is risky, but I couldn't get it to work without)
        else:
            subprocess.run([cmd] + args.split(), stdout=stdout, stderr=stderr)
    except subprocess.CalledProcessError:
        if not lenient:
            print("ERROR!! Cannot execute '%s'\n" % cmdline)
            exit()
        else:
            print("WARNING! Cannot execute '%s'\n" % cmdline)


def beamtilt(star_file=None, outdir='.', input_stack=None, coords_dir=None, box_size=None,
             micrograph_dir=None, particle_dir=None, verbosity=0, debug=False, options=None, usempi=False):
    """
    FUNCTION DESCRIPTION.

    Arguments:
        outdir : Output directory
        options : (list) Command-line options
        verbosity : How much information to write to screen
        usempi : (boolean) Whether to use Message Passing Interface
    """

    # Set up MPI
    mpioptions = setup_mpi(usempi)
    is_main = mpioptions['is_main']

    # Set output directory and log file name
    if is_main:
        prepare_outdir(outdir, main=is_main)
        log, _ = prepare_log(outdir)
    else:
        log = None

    quick_barrier()
    sp_global_def.write_command(
        outdir)  # assumes mpi_comm_rank=0, which isn't necessarily main process, which is random
    quick_barrier()

    do_continue = True
    if star_file:
        if is_main and verbosity >= 1:
            print_log_msg("Input STAR file: %s" % star_file, log)

    elif input_stack:
        # If not STAR file, need: tiltx, tilty, kv, cs, apix
        missing_items = []

        if not options.kv:
            do_continue = False
            missing_items.append("  --kv : Accelerating voltage, kV")

        # Some values can be 0, so "if not value" has to be "if value == None")
        if options.cs == None:
            do_continue = False
            missing_items.append("  --cs : Spherical aberration constant, mm")

        if not options.apix:
            do_continue = False
            missing_items.append("  --apix : Pixel size, Angstroms")

        if options.tiltx == None:
            do_continue = False
            missing_items.append("  --tiltx : Beam tilt along x, milliradians")

        if options.tilty == None:
            do_continue = False
            missing_items.append("  --tilty : Beam tilt along y, milliradians")

        if mpioptions['size'] != 1:
            do_continue = False
            if is_main:
                print_log_msg("ERROR! Parallelization not yet implemented without STAR file")

        if len(missing_items) > 0:
            msg = "ERROR! Without a STAR file, the following parameters must be supplied on the command line:"
            if is_main: print_log_msg(msg, log)

            if is_main:
                for item in missing_items: print_log_msg(item, log)

        else:
            if verbosity >= 1 and is_main:
                print_log_msg("Input stack: %s" % input_stack, log)
                print_log_msg("Accelerating voltage: %s kV" % options.kv, log)
                print_log_msg("Spherical aberration constant: %s mm" % options.cs, log)
                print_log_msg("Pixel size: %s Angstroms" % options.apix, log)
                print_log_msg("Beam tilt along x: %s milliradians" % options.tiltx, log)
                print_log_msg("Beam tilt along y: %s milliradians" % options.tilty, log)

    # Neither input stack nor STAR file
    else:
        do_continue = False
        if is_main:
            print_log_msg("ERROR!! Either STAR file or input stack must be supplied!", log)

    if not do_continue:
        if is_main: print()
        safe_exit()

    if is_main and verbosity >= 1:
        print_log_msg("Output directory: %s" % outdir, log)

        print_log_msg("Verbosity level set to %s (max %s)" % (verbosity, MAX_VERBOSE), log)
        print_log_msg("Debugging set to '%s'\n" % debug, log)

    # Read STAR file
    if star_file:
        star_obj = StarFile(star_file)
        star_data = star_obj.imported_content[""]  # string after "data_" tag, can't be None, can be ""
        star_parts = len(star_data)

        if not input_stack:
            if not coords_dir: coords_dir = os.path.join(outdir, EXTRACT_DIR)
            if verbosity >= 1 and is_main:
                msg = "Input stack not provided, will write particle coordinates to '%s/'" % coords_dir
                print_log_msg(msg, log)

        else:
            if is_main:
                if verbosity >= 1:
                    msg = "Input stack: %s" % input_stack
                    if debug: msg = "DEBUG " + msg
                    print_log_msg(msg, log)

                if not debug:
                    stack_parts = EMAN2.EMUtil.get_image_count(input_stack)

                    if star_parts != stack_parts:
                        print(
                            "ERROR!! Number of particles in '%s' and '%s' disagree (%s vs. %s)" % (
                                input_stack, star_file, stack_parts, star_parts
                            ))
                        exit()

        # If overriding STAR file values
        if verbosity >= 1 and is_main:
            if options.kv: print_log_msg("Accelerating voltage: %s kV" % options.kv, log)
            if options.cs: print_log_msg("Spherical aberration constant: %s mm" % options.cs, log)
            if options.apix: print_log_msg("Pixel size: %s Angstroms" % options.apix, log)
            if options.tiltx: print_log_msg("Beam tilt along x: %s milliradians" % options.tiltx, log)
            if options.tilty: print_log_msg("Beam tilt along y: %s milliradians" % options.tilty, log)
    else:
        star_data = None

    # Wait
    quick_barrier()

    # Extract particles if necessary (TODO: Move to function)
    if not input_stack:
        if is_main:
            if options.box_size == None:
                print("ERROR!! Box size required to extract coordinates! (flag '--box_size')")
                exit()

            cmd = 'sp_relion2sphire.py'
            args = '%s %s --box_size=%s --do_not_create_stack' % (
                star_file, os.path.join(outdir, EXTRACT_DIR), options.box_size
            )

            if verbosity >= 2:
                verbose = True
                stdout = None
            else:
                verbose = False
                stdout = subprocess.DEVNULL

            system_call_23(cmd, args, verbose=verbose, log=log, stdout=stdout)

            # Let user extract particles. Parallel process will be faster and shouldn't be spawned here.

            # Read STAR file data
            mic_path = star_data.loc[0, '_rlnMicrographName']
            pxsz = star_data.loc[0, '_rlnDetectorPixelSize']

            # Coordinates files from sp_relion2sphire
            box_template = os.path.join(coords_dir, os.path.dirname(mic_path), 'Coordinates/*.box')
            num_box_files = len(glob.glob(box_template))
            assert num_box_files > 0, "ERROR!! Couldn't find sp_relion2sphire output '%s'" % box_template

            # If user did not specify micrograph path, then get it from STAR file
            if micrograph_dir == None:
                micrograph_dir = os.path.dirname(mic_path)

            # Get micrograph extension
            mic_ext = os.path.splitext(os.path.basename(mic_path))[1]
            mic_template = os.path.join(micrograph_dir, "*%s" % mic_ext)
            if len(glob.glob(mic_template)) == 0: micrograph_dir = "<micrograph_directory>"

            # Use generic name for particle-extraction directory if not specified
            if particle_dir == None:
                particle_dir = "<output_particle_directory>"

            print()
            cmdline = "Particles ready to be extracted. "
            cmdline += "\nRun particle-extraction from the GUI, or from the command line using:"
            cmdline += "\n   sp_window.py '%s/*%s' '%s' %s %s --box_size=%s" % (
                micrograph_dir, mic_ext, box_template, pxsz, particle_dir, options.box_size
            )

            # Sanity checks
            if len(glob.glob(mic_template)) == 0:
                cmdline += "\nReplace '%s' above with the desired input directory name." % micrograph_dir
            if options.partdir == None:
                cmdline += "\nReplace '%s' above with the desired output directory name." % particle_dir

            print_log_msg(cmdline, log)

        print()
        safe_exit()
    # End extract IF-THEN

    quick_barrier()

    untilt_particles(
        input_stack,
        outdir,
        os.path.join(outdir, UNTILT_DIR),
        COMBINEDBDB,
        star_data=star_data,
        tilt_x=options.tiltx,
        tilt_y=options.tilty,
        voltage=options.kv,
        Cs=options.cs,
        angpix=options.apix,
        verbosity=verbosity,
        log=log,
        debug=debug,
        mpioptions=mpioptions
    )

    if is_main and verbosity > 0:
        print()
        print_log_msg("Done!\n", log)


def untilt_particles(input_stack, outdir, stack_dir, combined_bdb_stem, star_data=None,
                     tilt_x=None, tilt_y=None, voltage=None, Cs=None, angpix=None,
                     verbosity=0, log=None, debug=False, mpioptions=None):
    is_main = mpioptions['is_main']
    usempi = mpioptions['use']
    mpi_id = mpioptions['myid']

    # Separate by micrograph
    if isinstance(star_data, pandas.core.frame.DataFrame):
        micrograph_dict = separate_by_mic(
            star_data,
            verbosity=verbosity,
            log=log,
            debug=debug,
            mpioptions=mpioptions
        )
        num_tot = len(star_data)
    else:
        micrograph_dict = {}
        num_tot = EMAN2.EMUtil.get_image_count(input_stack)
        micrograph_dict[input_stack] = list(range(num_tot))

    quick_barrier()

    if is_main:
        # Clean up pre-existing output
        if not os.path.exists(stack_dir): os.makedirs(stack_dir)

    quick_barrier()

    # Don't calculate voltage every time, don't assume voltage is constant
    prevVolts = -1

    if usempi:
        mpi_start, mpi_end = sp_applications.MPI_start_end(
            len(micrograph_dict),
            mpi.mpi_comm_size(mpi.MPI_COMM_WORLD),
            mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        )
    else:
        mpi_start = 0
        mpi_end = len(micrograph_dict)

    # Make a subset of the dictionary
    mpi_dict = collections.OrderedDict(itertools.islice(micrograph_dict.items(), mpi_start, mpi_end))

    if is_main and verbosity >= 1:
        if not debug:
            msg = "Applying beam-tilt correction to %s images in '%s'" % (num_tot, input_stack)
            print_log_msg(msg, log)
        else:
            print_log_msg("DEBUG: Not applying corrections", log)

    if verbosity >= 1:
        if usempi:
            if is_main:
                msg = "Correcting %s micrograph stacks in chunks of %s" % (len(micrograph_dict), len(mpi_dict))
                print_log_msg(msg, log)
        else:
            print_log_msg("Correcting %s micrograph stacks in series" % (len(micrograph_dict)), log)

    quick_barrier()

    disable_mic = verbosity < 2 or verbosity > 2 or debug or not is_main

    # Loop through micrographs
    for mic_idx, micname in tqdm.tqdm(
            enumerate(mpi_dict),
            unit='mic',
            disable=disable_mic,
            file=sys.stdout
    ):

        # Set up filenames
        mic_stack = os.path.join(stack_dir, micname)
        if os.path.exists(mic_stack): os.remove(mic_stack)

        # Set up BDB name
        mic_stem = os.path.splitext(micname)[0]
        mic_bdb_obj = BdbNames(mic_stem, outdir)
        if os.path.exists(mic_bdb_obj.bdb_path): os.remove(mic_bdb_obj.bdb_path)  # or else will append
        mic_bdb_dict = EMAN2db.db_open_dict(mic_bdb_obj.bdb_name)

        if is_main and verbosity >= 3:
            print_log_msg("Working on stack %s (%s/%s)" % (micname, mic_idx + 1, len(mpi_dict)), log)

        disable_part = verbosity < 3 or debug or not is_main

        # Loop through particles
        for part_idx in tqdm.tqdm(
                range(len(micrograph_dict[micname])),
                unit='img',
                disable=disable_part,
                file=sys.stdout):

            glonum = micrograph_dict[micname][part_idx]

            # Some values can be 0, so "if not value" has to be "if value == None")
            if not voltage: voltage = star_data.loc[glonum, '_rlnVoltage']
            if Cs == None: Cs = star_data.loc[glonum, '_rlnSphericalAberration']
            if not angpix: angpix = star_data.loc[glonum, '_rlnDetectorPixelSize']
            if tilt_x == None: tilt_x = star_data.loc[glonum, '_rlnBeamTiltX']
            if tilt_y == None: tilt_y = star_data.loc[glonum, '_rlnBeamTiltY']

            if voltage != prevVolts:
                wavelength = 1e10 * getElectronWavelength(voltage * 1000)  # convert m to Angstroms, kV to V

            if verbosity >= 6:
                msg = "glonum %s, micname %s, tilt_x %s, tilt_y %s, wavelength %.4f, Cs %s, angpix %s" % (
                    glonum, micname, tilt_x, tilt_y, wavelength, Cs, angpix
                )
                print_log_msg(msg, log)

            if not debug:
                im0 = sp_utilities.get_im(input_stack, glonum)

                # Apply beam-tilt correction
                corr_img, _ = apply_tilt(
                    im0,
                    tilt_x,
                    tilt_y,
                    wavelength,
                    Cs,
                    angpix,
                    log=log,
                    verbose=verbosity >= 7
                )

                corr_img.write_image(mic_stack, part_idx)

                # Copy header
                header_dict = im0.get_attr_dict()
                header_dict["ptcl_source_coord_id"] = part_idx
                header_dict["data_path"] = sp_utilities.makerelpath(mic_bdb_obj.eman2db_dir, mic_stack)
                mic_bdb_dict[part_idx] = header_dict
            else:
                pathlib.Path(mic_bdb_obj.bdb_path).touch()  # placeholder
        # End particle loop

        # Close BDB
        EMAN2db.db_close_dict(mic_bdb_obj.bdb_name)
    # End micrograph loop

    quick_barrier()

    if is_main:
        if verbosity >= 2 or debug: print()

        # Merge BDBs
        merge_bdbs(
            outdir,
            outdir,
            combined_bdb_stem,
            verbosity=verbosity,
            log=log,
            debug=debug
        )

    quick_barrier()


def separate_by_mic(star_data, verbosity=0, log=None, debug=False, mpioptions=None):
    """

    Returns:
        micrograph_dict : Dictionary containing:
            key : stack name, i.e., for each micrograph
            list of particles
    """

    # Only the main process should print anything
    if mpioptions == None:
        verbosity = 0
    else:
        if not mpioptions['is_main']: verbosity = 0

    if verbosity >= 1: print_log_msg("Separating by micrograph...", log)

    micrograph_dict = collections.OrderedDict()

    # Loop through particles (fast, don't need progress bar)
    for idx in range(len(star_data)):
        micname = os.path.basename(star_data.loc[idx, '_rlnImageName'])

        # Create an empty list for each new micrograph
        if micname not in micrograph_dict.keys(): micrograph_dict[micname] = []

        # Append particle index to list
        micrograph_dict[micname].append(idx)

    if verbosity >= 5:
        for micname in micrograph_dict.keys():
            print_log_msg("Stack %s: Found %s particles" % (
                micname, len(micrograph_dict[micname])
            ), log)

    if verbosity >= 1: print_log_msg("Separated %s particles into %s micrographs\n" % (
        len(star_data), len(micrograph_dict)
    ), log)

    return micrograph_dict


def getElectronWavelength(voltage):
    # voltage in Volts, length unit in meters
    # Adapted from Anchi Cheng's fftfun.py

    h = 6.6e-34
    m = 9.1e-31
    charge = 1.6e-19
    c = 3e8
    wavelength = h / math.sqrt(2 * m * charge * voltage)
    relativistic_correction = 1 / math.sqrt(1 + voltage * charge / (2 * m * c * c))

    return wavelength * relativistic_correction


def apply_tilt(img, tilt_x, tilt_y, wavelength, Cs, angpix, mode='scipy', log=None, verbose=False, debug=None):
    assert img['nx'] == img['ny'], "ERROR!! Assuming square images! %s x %s" % (img['nx'], img['ny'])
    # (I don't want to deal with unequal resolution in x & y. Maybe we could pad the image to a square.)

    # Convert to numpy
    imgarray = img.get_2dview()

    if mode == 'scipy':
        imgfft = scipy.fft.rfft2(imgarray)
        ####print('imgfft', type(imgfft))
        imgfft = np.transpose(imgfft)
        ####print('imgfft', imgfft.shape)
    else:
        imgfft = np.fft.fft2(imgarray)

    xdim = imgfft.shape[0]
    ydim = imgfft.shape[1]
    ####printvars(['xdim','ydim'], True)

    boxsize = angpix * ydim
    centerx = xdim // 2 + 1
    centery = ydim // 2 + 1

    # Convert from milliradians to degrees
    factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize)

    if verbose: print("xc yc jp ip   magnitd  old_phase      delta   new_real  new_phase")
    if debug == 'spider': print("jp ip   oldreal    oldimag      delta    newreal    oldimag")

    if mode == 'scipy':
        # Loop through y
        for ycoord in range(ydim):
            if ycoord < centery:
                ip = ycoord
            else:
                ip = ycoord - ydim

            # Loop through x
            for xcoord in range(xdim):
                jp = xcoord

                # Calculate phase shift in degrees
                delta_phase = factor * (ip * ip + jp * jp) * (ip * tilt_y + jp * tilt_x)

                realval = imgfft[xcoord, ycoord].real
                imagval = imgfft[xcoord, ycoord].imag
                mag = np.sqrt((realval * realval) + (imagval * imagval))
                old_phase = np.arctan2(imagval, realval)

                # Convert phase shift from degrees to radians
                new_phase = old_phase + np.deg2rad(delta_phase)

                newrealval = mag * np.cos(new_phase)
                newimagval = mag * np.sin(new_phase)

                imgfft[xcoord, ycoord] = newrealval + (1j * newimagval)

                if verbose:
                    print(
                        f"{xcoord}, {ycoord}, {jp}, {ip}, {mag: 9.4f}, {np.rad2deg(old_phase): 9.4f}, {delta_phase: 9.4f}, {newrealval: 9.4f}, {newimagval: 9.4f}")

                if debug == 'spider':
                    print(
                        f'{jp}, {ip}, {realval: 9.4f}, {imagval: 9.4f}, {delta_phase: 9.4f}, {newrealval: 9.4f}, {newimagval: 9.4f}')

        imgfft = np.transpose(imgfft)
        real_array = scipy.fft.irfft2(imgfft)

    else:
        # Loop through y
        for ycoord in range(ydim):
            if ycoord < centery:
                ip = ycoord
            else:
                ip = ycoord - ydim

            # Loop through x
            for xcoord in range(centerx):
                jp = xcoord

                # Calculate phase shift in degrees
                delta_phase = factor * (ip * ip + jp * jp) * (ip * tilt_y + jp * tilt_x)

                realval = imgfft[ycoord, xcoord].real
                imagval = imgfft[ycoord, xcoord].imag
                mag = np.sqrt((realval * realval) + (imagval * imagval))
                old_phase = np.arctan2(imagval, realval)

                # Convert phase shift from degrees to radians
                new_phase = old_phase + np.deg2rad(delta_phase)

                newrealval = mag * np.cos(new_phase)
                newimagval = mag * np.sin(new_phase)

                imgfft[ycoord, xcoord] = newrealval + (1j * newimagval)

                if verbose:
                    print(
                        f"{xcoord}, {ycoord}, {jp}, {ip}, {mag: 9.4f}, {np.rad2deg(old_phase): 9.4f}, {delta_phase: 9.4f}, {newrealval: 9.4f}, {newimagval: 9.4f}")

                if debug == 'spider':
                    print(
                        f'{jp}, {ip}, {realval: 9.4f}, {imagval: 9.4f}, {delta_phase: 9.4f}, {newrealval: 9.4f}, {newimagval: 9.4f}')

        real_array = np.fft.ifft2(imgfft)
    # End scipy if-then

    eman_img = sp_utilities.numpy2em_python(real_array)  # .real)  # EMNumPy.numpy2em(real_array.real)

    return eman_img, np.transpose(imgfft)


def merge_bdbs(input_bdb_dir, output_bdb_dir, output_bdb_stem, verbosity=0, log=None, debug=False):
    """
    Merges a series of BDBs.
    Adapted from sp_tilt_import.py

    Modified 2020-11-26

    Arguments:
        input_bdb_dir : directory for input BDBs
        output_bdb_dir : target directory for merged BDB (not including EMAN2DB subdirectory)
        output_bdb_stem : stem of output BDB (that is, without preceding "bdb:" or trailing ".bdb")
        verbosity : how much information to write to screen
        log : instance of Logger class
        debug : (boolean) for debugging
    """

    # Output BDB
    output_bdb_obj = BdbNames(output_bdb_stem, output_bdb_dir)
    if os.path.exists(output_bdb_obj.bdb_path): os.remove(
        output_bdb_obj.bdb_path)  # will otherwise merge with pre-existing file

    if verbosity >= 1:
        print()
        print_log_msg("Combining BDBs from '%s' into '%s'" % (input_bdb_dir, output_bdb_obj.bdb_name), log)

    input_eman_bdb_pattern = os.path.join(input_bdb_dir, 'EMAN2DB', "*.bdb")
    input_bdb_list = sorted(glob.glob(input_eman_bdb_pattern))

    # Remove files that shouldn't be merged
    bad_list = ([s for s in input_bdb_list if '00image_counts.bdb' in s])
    for bad_file in bad_list:
        if verbosity >= 3 or debug: print_log_msg('Excluding %s from merging step' % bad_file, log)
        input_bdb_list.remove(bad_file)

    if debug:
        if verbosity > 0:
            for input_bdb_file in input_bdb_list:
                msg = "DEBUG: Will merge %s into %s" % (input_bdb_file, output_bdb_obj.bdb_name)
                print_log_msg(msg, log)
        return

    # We'll need this variable when updating the data_path in the header
    out_eman2db_dir = os.path.abspath(output_bdb_obj.eman2db_dir)

    # Open new database
    new_bdb_dict = EMAN2db.db_open_dict(output_bdb_obj.bdb_name)

    img_counter = 0

    # Turn off progress bar if we're already printing information every iteration
    disableTF = verbosity < 2 or verbosity > 3

    # Loop through classes
    for input_bdb_file in tqdm.tqdm(input_bdb_list, unit=' BDB file', disable=disableTF, file=sys.stdout):
        input_bdb_stem = os.path.splitext(os.path.basename(input_bdb_file))[0]
        input_bdb_obj = BdbNames(input_bdb_stem, input_bdb_dir)
        input_bdb_name = input_bdb_obj.bdb_name

        input_bdb_dict = EMAN2db.db_open_dict(input_bdb_name, ro=True)

        # Get number of particles
        bdb_numparts = len(input_bdb_dict)
        if verbosity >= 4:
            print_log_msg('Opened %s images from %s' % (bdb_numparts, input_bdb_name), log)

        # Loop through particles
        for partnum in range(bdb_numparts):
            # Read image header
            img_header = input_bdb_dict.get_header(partnum)

            # Update relative path between new BDB and image stack
            rel_path = update_data_path(
                img_header['data_path'],
                os.path.join(input_bdb_dir, 'EMAN2DB'),
                out_eman2db_dir)
            img_header['data_path'] = rel_path

            # (When are the following used?)
            img_header['data_source'] = input_bdb_name
            img_header['data_n'] = partnum

            # Write in new database
            new_bdb_dict[img_counter] = img_header
            img_counter += 1
        # End particle loop

        # Clean up
        EMAN2db.db_close_dict(input_bdb_name)
        del input_bdb_obj
    # End class loop

    # Close new database
    EMAN2db.db_close_dict(output_bdb_obj.bdb_name)
    numoutimgs = EMAN2.EMUtil.get_image_count(output_bdb_obj.bdb_name)
    assert numoutimgs == img_counter, "Uh oh!! merge_bdbs: %s != %s" % (numoutimgs, img_counter)
    if verbosity >= 1:
        print_log_msg('Wrote %s images to %s' % (numoutimgs, output_bdb_obj.bdb_name), log)


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
    bdb2curr_path = os.path.join(in_eman2db_dir, curr_path)

    # Get absolute directory
    abs_img_path = os.path.abspath(bdb2curr_path)
    # (abspath removes '..'s from the path.)

    # Sanity check
    if not os.path.exists(abs_img_path):
        print("ERROR!! update_data_path: %s does not exist!" % abs_img_path)
        printvars(['in_eman2db_dir', 'curr_path', 'out_eman2db_dir', 'abs_img_path'], quitTF=True)

    # Get relative path between source (BDB) and target (image stack)
    rel_path = sp_utilities.makerelpath(out_eman2db_dir, abs_img_path)

    return rel_path


def parse_command_line():
    """
    Parse the command line.  Adapted from sxmask.py

    Arguments:
    None:

    Returns:
    Parsed arguments object
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        usage=USAGE,
        epilog=MODIFIED
    )

    parser.add_argument(
        'outdir',
        type=str,
        help='Output directory')

    parser.add_argument(
        '--star_file',
        type=str,
        help='Input STAR file')

    parser.add_argument(
        '--input_stack',
        type=str,
        default=None,
        help='Input stack')

    parser.add_argument('--verbosity', "-v",
                        type=int,
                        default=3,
                        help='Increase verbosity')

    """
    Verbosity levels:
        0: None
        1: Basic
        2: Progress bar for micrographs
        2: Display sp_relion2sphire output
        3: Particles per micrograph
        3: Progress bar for particles within a micrograph
        3: BDB files excluded from merging
        4: Merged individual BDBs
        5: Particles separates by micrograph
        6: Stats for each particle
        7: For each pixel (!)
    """

    parser.add_argument('--debug',
                        action="store_true",
                        help='Dry run, for debugging purposes')

    group_extract = parser.add_argument_group(
        title='Extraction options, for use in sp_window call',
        description="Options if input stack doesn't exist.")

    group_extract.add_argument(
        '--micdir',
        type=str,
        default=None,
        help='Input directory where micrographs are located')

    group_extract.add_argument(
        '--box_size',
        type=int,
        default=None,
        help='Dimension for extracted particles')

    group_extract.add_argument(
        '--coordsdir',
        type=str,
        default=None,
        help='Output directory where particle coordinates will be written')

    group_extract.add_argument(
        '--partdir',
        type=str,
        default=None,
        help='Output directory where particles will be written')

    group_nostar = parser.add_argument_group(
        title='Fixed parameters independent of STAR file.',
        description="Parameters if not using, or overriding, STAR file")

    group_nostar.add_argument(
        '--kv',
        type=float,
        default=None,
        help='Accelerating voltage, kV')

    group_nostar.add_argument(
        '--cs',
        type=float,
        default=None,
        help='Spherical aberration constant, mm')

    group_nostar.add_argument(
        '--apix',
        type=float,
        default=None,
        help='Pixel size, Angstroms')

    group_nostar.add_argument(
        '--tiltx',
        type=float,
        default=None,
        help='Beam tilt along x, milliradians')

    group_nostar.add_argument(
        '--tilty',
        type=float,
        default=None,
        help='Beam tilt along y, milliradians')

    return parser.parse_args()


def main():
    options = parse_command_line()

    ##print args, options  # (Everything is in options.)
    # print(options)
    ##print('LOGFILE',sp_global_def.LOGFILE)
    # exit()

    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
        mpi.mpi_init(0, [])  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp("Start")

    beamtilt(
        star_file=options.star_file,
        outdir=options.outdir,
        input_stack=options.input_stack,
        coords_dir=options.coordsdir,
        box_size=options.box_size,
        micrograph_dir=options.micdir,
        particle_dir=options.partdir,
        verbosity=options.verbosity,
        debug=options.debug,
        options=options,
        usempi=RUNNING_UNDER_MPI
    )

    sp_global_def.print_timestamp("Finish")
    if RUNNING_UNDER_MPI:
        mpi.mpi_finalize()


if __name__ == "__main__":
    main()
