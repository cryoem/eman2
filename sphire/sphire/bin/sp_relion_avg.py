#!/usr/bin/env python
from __future__ import print_function
from past.utils import old_div

from sphire.libpy import sp_global_def 
from sphire.libpy import sp_logger
from sphire.libpy import sp_fundamentals 
from sphire.libpy import sp_filter 
from sphire.libpy import sp_morphology
from sphire.libpy import sp_utilities
from sphire.libpy import sp_alignment
import mpi
import os
import datetime
import argparse
import inspect
import six
import re
import math
import tqdm
import EMAN2
import EMAN2_cppwrap
import sys
import copy
import pandas
import subprocess
import sys

TIMESTAMP_LENGTH= 23

CLASSDOCDIR=    "ClassDocs"    # Class text file directory
CLASSDOCPREFIX= "sel_class_"   # Class-membership files
CLASSAVGSTACK=  "stkavgs.hdf"  # Class average stack

USAGE = """ 
PURPOSE:
Computes averages using alignment parameters from STAR files.

%s <input_star_file> <output_directory>
Required command-line parameters:
    1. Input STAR file
    2. Output directory
Outputs:
    stkavgs.hdf : Stack of averages
    ClassDocs/sel_class_* : Selection files for each class

General options:
%s <input_star_file> <output_directory> --stackdir=<particle_stack_directory> --verbosity=<verbosity>
Parameters:
    --stackdir : Directory where input images are located, relative to path in STAR file
    --ctf : CTF-correction options: none (default), flip, wiener, punish
    --verbosity : Increase verbosity (0..5)
Advanced parameters:
    --negate : Negates signs for shifts read from STAR file
    --unaligned : STAR file before alignment (for debugging)

""" % ((__file__,)*2)

MODIFIED="Modified 2021-02-22"

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
            skiprows=line_dict['content'][0]-1,
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
            skiprows=line_dict['content'][0]-1,
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
    
    Argument:
        usempi : (boolean) Whether to use Message Passing Interface
    """
    
    if usempi:
        number_of_proc = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
        myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
        
        #  choose a random node as main
        main_node = 0
        if(myid == 0): main_node = random.randint(0,number_of_proc-1)
        main_node = mpi.mpi_bcast(main_node, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)	
    else:
        myid = main_node = 0
        # (IF-THEN statements will be simpler if I pretend non-MPI runs are the main node.)
    
    is_main = myid==int(main_node)  # boolean: whether current process is the main one
    mpidict = {'use':usempi, 'myid':myid, 'main_node':main_node, 'is_main':is_main}
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
    
    logname = "log_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") +  ".txt"
    logname = os.path.join(outdir, logname)
    
    # May be using old version of logger.py
    try:
        if verbose:
            log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), base_logger2=sp_logger.BaseLogger_Print(), file_name=logname)
            verbose = False  # logger output will be echoed to screen
        else:
            log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files(), file_name=logname)
    except TypeError:
        print("WARNING: Using old logger.py library")
        log = sp_logger.Logger(base_logger=sp_logger.BaseLogger_Files())#, file_name=logname)
        logname = 'log.txt'
        
    if main:
        print("Writing log file to %s" % logname)
    
        progbase = os.path.basename(__file__).split('.')[0].upper()
        length = len(progbase) + 4
        
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
    
    frame = inspect.currentframe()
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

def system_call_23(cmd, args, stdout=None, stderr=None, log=None, verbose=False):
    """
    NOT WORKING YET
    """
    
    # Check Python version
    python_version= sys.version.split()[0]
    
    # Keep up to one decimal
    version_float= float('.'.join(python_version.split('.', 2)[:2]))
    
    # Build command line
    cmdline= "%s %s" % (cmd, args)
    
    if verbose : print_log_msg(cmdline, log)
    
    if version_float < 3.5:
        subprocess.check_call([cmd, args], stdout=stdout, stderr=stderr, shell=True)
        print("UNTESTED!! Confirm success!")
        exit()
    else:
        subprocess.run(cmd, args, stdout=stdout, stderr=stderr)
        print("UNTESTED!! Confirm success!")
        exit()

def relion_avg(input_star_file, outdir='.', options=None, verbosity=0, usempi=False):
    """
    FUNCTION DESCRIPTION.
    
    Arguments:
        outdir : Output directory
        options : (list) Command-line options
        verbosity : How much information to write to screen
        usempi : (boolean) Whether to use Message Passing Interface
    """
    
    if hasattr(options, 'unaligned'):
        unaligned_star= options.unaligned
    else:
        unaligned_star= None
    
    if hasattr(options, 'stackdir'):
        stack_dir= options.stackdir
    else:
        stack_dir= None
    
    #Set up MPI
    mpioptions = setup_mpi(usempi)
    
    # Set output directory and log file name
    prepare_outdir(outdir, main=mpioptions['is_main'])
    sp_global_def.write_command(outdir)
    log, _ = prepare_log(outdir, main=mpioptions['is_main'])
    
    # Wait
    if usempi: mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    
    if verbosity>= 1:
        print_log_msg('Input STAR file: %s' % input_star_file, log)
        print_log_msg('Output directory: %s' % outdir, log)
        
        if unaligned_star:
            print_log_msg('Will compare alignment parameters to: %s' % unaligned_star, log)
        
        if stack_dir:
            print_log_msg('Will assume particle stacks to be in %s/' % stack_dir, log)
        
        print_log_msg("CTF-correction mode set to '%s'" % options.ctf, log)
        
        if options.negate:
            print_log_msg('Will negate shifts in STAR file', log)
        
        print_log_msg('Verbosity level (max 5): %s\n' % verbosity, log)
        
    if verbosity>= 1:
        print_log_msg('Reading STAR file: %s' % input_star_file, log)
    
    # Read STAR file
    final_star_obj= StarFile(input_star_file)
    final_data= final_star_obj.imported_content['particles']
    num_parts= len(final_data)
    apix= float(final_star_obj.imported_content['optics']['_rlnImagePixelSize'])
    final_dict= final_data[0: num_parts]
    class_list= final_dict['_rlnClassNumber']
    
    if verbosity>= 1:
        print_log_msg('Finished reading STAR file with %s entries\n' % num_parts, log)
    
    # Compare old angles to new, make sure only psi differs (TODO: Move to function)
    if unaligned_star: 
        # Read data into lists
        init_star_obj= StarFile(unaligned_star)
        init_data= init_star_obj.imported_content['particles']
        
        # Sanity check
        assert len(init_data) == num_parts, "ERROR!! Initial and final STAR files have different legnths: %s vs. %s" % (len(init_data), num_parts)
        
        if verbosity>=1 : print_log_msg('Comparing initial to final angles', log)
        
        # Loop through particles
        for idx in range(num_parts):
            # Read final data
            curr_rot_final= final_data.loc[idx, '_rlnAngleRot']
            curr_tilt_final= final_data.loc[idx, '_rlnAngleTilt']
            
            # Read initial data
            curr_rot_init= init_data.loc[idx, '_rlnAngleRot']
            curr_tilt_init= init_data.loc[idx, '_rlnAngleTilt']
            
            # Sanity check
            if curr_rot_final != curr_rot_init or curr_tilt_final != curr_tilt_init:
                print("ERROR!! Unknown state!")
                printvars([
                    'idx', 
                    'curr_rot_final', 
                    'curr_rot_init', 
                    'curr_tilt_final', 
                    'curr_tilt_init'
                    ], quitTF=True)
        # End particle loop    
        
        if verbosity>=1 : print_log_msg('Initial and final angles OK...\n', log)
    # End unaligned IF-THEN
    
    # Get class information (TODO: Move to function)
    class_doc_dir= os.path.join(outdir, CLASSDOCDIR)
    if not os.path.isdir(class_doc_dir) : os.makedirs(class_doc_dir)
    
    if verbosity>=1 : print_log_msg('Checking classes', log)
    
    class_dict= {}
    
    # Loop through particles
    for idx in range(num_parts):
        curr_class= class_list[idx]
        
        if not (curr_class in class_dict):
            class_dict[curr_class]= [idx]
        else:
            class_dict[curr_class].append(idx)
    # End particle loop
    
    populated_classes= list(class_dict.keys() )
    populated_classes.sort()
    num_classes= len(populated_classes)
    num_digits= math.ceil( math.log10(num_classes) )
    digit_pattern= '{{0:0{0}d}}'.format(num_digits)
    
    tot_class_parts= 0
    sphire2relion_lut= {}
    
    # Print class population
    for class_num in range(num_classes):
        relion_class= populated_classes[class_num]
        tot_class_parts+= len(class_dict[relion_class])
        
        # Map RELION class to SPHIRE class
        sphire2relion_lut[class_num]= relion_class
        
        # Write members to file
        class_doc_file= os.path.join(class_doc_dir, CLASSDOCPREFIX) + digit_pattern.format(class_num) + '.txt'
        sp_utilities.write_text_row(class_dict[relion_class], class_doc_file)
        
        if verbosity>=4 : 
            print_log_msg("%s:  Wrote %s members from RELION class %s to %s" % (
                class_num, len(class_dict[relion_class]), relion_class, class_doc_file
                ), log)
    # End class loop
    
    if verbosity>=1: 
        print_log_msg('Assigned %s particles to %s classes and averaging\n' % (tot_class_parts, num_classes), log)
    
    # Apply shifts (TODO: move to function)
    
    if verbosity>=1 : print_log_msg('Applying shifts', log)
    
    disable_classes= verbosity>1 or verbosity<1
    
    # Loop through classes
    for class_num in tqdm.tqdm(range(num_classes), unit='class', disable=disable_classes, file=sys.stdout):
        # Get RELION class number
        relion_class= sphire2relion_lut[class_num]
        
        aligned_list= []
        
        disable_parts= verbosity>4 or verbosity<3
        
        if verbosity>=2:
            msg= "Class %s: Relion class %s, %s particles" % (class_num, relion_class, len(class_dict[relion_class]) )
            print_log_msg(msg, log)

        # Loop through particles
        for part_idx in tqdm.tqdm(range( len(class_dict[relion_class]) ), unit='img', disable=disable_parts, file=sys.stdout):
            part_num= class_dict[relion_class][part_idx]
            
            # Read shifts, rotation, etc.
            curr_psi=      final_data.loc[part_num, '_rlnAnglePsi']
            curr_xshift=   final_data.loc[part_num, '_rlnOriginXAngst']/apix
            curr_yshift=   final_data.loc[part_num, '_rlnOriginYAngst']/apix
            curr_part2mic= final_data.loc[part_num, '_rlnImageName']
            curr_class=    class_list[part_num]
            
            # Separate slice number from stack name
            slice_num, stack_path= curr_part2mic.split("@")
            slice_num= int(slice_num)
            
            if stack_dir:
                ####final_stack= os.path.join(stack_dir, os.path.basename(stack_path) )
                final_stack= os.path.join(stack_dir, stack_path)
            else:
                final_stack= stack_path
            
            # Initialize stack if necessary
            if part_idx == 0 and class_num == 0:
                # Read dimensions 
                box_size= sp_utilities.get_im(final_stack, slice_num-1)['nx']
                
                # Initialize EMData object to write aligned particles
                aligned_stack_obj= EMAN2.EMData(box_size, box_size, num_classes)  # (TODO: Don't treat as volume)
                
            if verbosity>=5:
                msg= "Class %s, particle %s, angle %s, shift (%s, %s), stack %s@%s" % (
                    curr_class, part_num, curr_psi, curr_xshift, curr_yshift, slice_num, final_stack)
                print_log_msg(msg, log)
                exit()
                
            sphire_cter_entry = {}

            relion_defocusU = float(final_data.loc[part_num, "_rlnDefocusU"])
            relion_defocusV = float(final_data.loc[part_num, "_rlnDefocusV"])
            relion_defocus_angle = float(final_data.loc[part_num, "_rlnDefocusAngle"])

            sphire_defocus= (relion_defocusU + relion_defocusV)/20000  # convert format from RELION to SPHIRE
            astig_amp= (-relion_defocusU + relion_defocusV)/10000  # convert format from RELION to SPHIRE
            
            astig_ang= 45.0 - relion_defocus_angle  # convert format from RELION to SPHIRE
            while astig_ang >= 180 : astig_ang-= 180
            while astig_ang < 0 : astig_ang+= 180

            const_ac= 100. * float(final_star_obj.imported_content['optics']["_rlnAmplitudeContrast"])
            phase_shift = float(final_data.loc[part_num, "_rlnPhaseShift"])
            ac_phase_shift = sp_morphology.ampcont2angle(const_ac)  # must pass amplitude contrast in [%]
            total_phase_shift = (phase_shift + ac_phase_shift)
            total_ac = sp_morphology.angle2ampcont(total_phase_shift)
            # (NOTE: Untested with phase shift)

            voltage= float(final_star_obj.imported_content['optics']['_rlnVoltage'])
            spher_aberr= float(final_star_obj.imported_content['optics']['_rlnSphericalAberration'])
            
            # Sanity check
            part_name= int(final_data.loc[part_num, "_rlnOriginalParticleName"])
            assert part_name == part_num, "Uh oh! Particle numbers don't match (%s vs %s)" % (part_name, part_num)

            ctf_obj= sp_utilities.generate_ctf([
                sphire_defocus, 
                spher_aberr, 
                voltage, 
                apix, 
                0.0,
                total_ac, 
                astig_amp, 
                astig_ang
                ])
            
            # Read image
            try:
                img_orig= sp_utilities.get_im(final_stack, slice_num-1)
            except RuntimeError:
                printvars([
                    'final_stack',
                    'slice_num'
                    ], quitTF=True)
            
            # Apply alignment parameters (default is not to negate shifts, see '--negate')
            if not options.negate:
                img_ali= sp_fundamentals.rot_shift2D(img_orig, -curr_psi,  curr_xshift,  curr_yshift, interpolation_method='quadratic')
            else:
                img_ali= sp_fundamentals.rot_shift2D(img_orig, -curr_psi, -curr_xshift, -curr_yshift, interpolation_method='quadratic')
            
            # CTF parameters
            img_ali.set_attr("ctf", ctf_obj)
            
            # Write to stack
            aligned_list.append(img_ali)
        # End particle loop
        
        # Compute average
        ####avg_img, _= avgvar_ctf(aligned_list, dopa=True, mode='n')  
        avg_img, _= avg_optional_ctf(aligned_list, ctf_method=options.ctf, do_align=False, progressbar=(not disable_parts) )  
        # 'dopa' pads, 'n' does not apply alignment
        
        avg_img.set_attr('relion_class', int(relion_class) )
        aligned_stack_obj.insert_clip(avg_img, (0, 0, class_num) )
    # End class loop
    
    avg_stack_file= os.path.join(outdir, CLASSAVGSTACK)
    aligned_stack_obj.write_image(avg_stack_file)
    
    if verbosity>=1 : 
        print_log_msg('Finished applying shifts', log)
        print_log_msg('Wrote averages to %s' % avg_stack_file , log)
    
def avg_optional_ctf(bdb_or_list, ctf_method=None, pws_docfile=None, do_align=True, progressbar=False):
    """
    
    Arguments:
        bdb_or_list : input image BDB stack or list of EMData objects
        ctf_method : method for CTF correctioni: None, 'flip', 'wiener'
        pws_doc_template : if phase-flipping, 1D profile of amplitudes before & after CTF-correction
        do_align : (boolean) whether to apply alignment to images in input stack
        progressbar : (boolean) whether to show progress bar
    
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
        avg_img, var_img= avgvar_ctf(bdb_or_list, dopa=True, mode=align_mode, progressbar=progressbar)
    
    elif ctf_method == 'flip':
        avg_img, var_img, pxsz= avgvar_flipctf(bdb_or_list, dopa=False, mode=align_mode, progressbar=progressbar)
        avg_noctf, var_noctf=   avgvar(bdb_or_list, mode=align_mode, progressbar=progressbar)
        
        # Compute 1D power spectra
        rops_ctf= sp_fundamentals.rops_table(avg_img)
        rops_noctf= sp_fundamentals.rops_table(avg_noctf)
        spfreq= [ x/(2.0 * pxsz * (len(rops_ctf)-1) ) for x in range(len(rops_ctf))]  # float(2) ensures float result
        pws_tuple_list= list( zip(spfreq, rops_ctf, rops_noctf) )  # list() should work in both python 2 & 3
        pws_list= map(list, pws_tuple_list)  # wrote_text_row won't like list of tuples
        
        # Write 1D power spectra
        if pws_docfile : sp_utilities.write_text_row(pws_list, pws_docfile)
    
    elif ctf_method == 'punish':
        # Apply trapped CTF to images
        if do_align: 
            print("ERROR!! There is currently no option for applying alignments with option '%s'" % ctf_method)
            exit()
        else:
            punished_img_list= punish_stack(bdb_or_list, progressbar=progressbar)
        avg_img, var_img= avgvar(bdb_or_list, mode=align_mode, progressbar=progressbar)
    
    elif ctf_method == None:
        avg_img, var_img= avgvar(bdb_or_list, mode=align_mode, progressbar=False)
    
    else:
        print("ERROR!! ctf_method %s not known" % ctf_method)
        exit()
    
    return avg_img, var_img

def avgvar(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, progressbar=False):
    '''
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
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
    inmem = True
    if type(data) == type(""):
        inmem = False

    img2D = True
    if inmem:
        img = data[0]
    else:
        img = sp_utilities.get_im(data,0)
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    if nz > 1:
        img2D = False

    if inmem:
        data_nima = len(data)
    else:
        data_nima = EMUtil.get_image_count(data)
    if i2 == 0: i2 = data_nima-1

    ave = sp_utilities.model_blank(nx,ny,nz)
    var = sp_utilities.model_blank(nx,ny,nz)
    nima = 0
    for i in tqdm.tqdm(range(i1, i2+1), unit='img', disable=(not progressbar), file=sys.stdout):
        if not(use_odd) and i%2 == 1:
            continue
        if not(use_even) and i%2 == 0:
            continue
        nima += 1
        if inmem:
            img = data[i]
        else:
            img = sp_utilities.get_im(data, i)
        
        #printvars('mode', quitTF=True)  #### DIAGNOSTIC
        
        if (mode == 'a'):
            if img2D:
                angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
                img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale, interp)
            else:
                phi, theta, psi, s3x, s3y, s3z, mirror, scale = sp_utilities.get_params3D(img)
                img = sp_fundamentals.rot_shift3D(img, phi, theta, psi, s3x, s3y, s3z, scale)
        EMAN2_cppwrap.Util.add_img(ave, img)
        EMAN2_cppwrap.Util.add_img2(var, img)

    EMAN2_cppwrap.Util.mul_scalar(ave, 1.0 /float(nima) )
    return ave, (var - ave*ave*nima)/(nima-1)

def avgvar_ctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0, dopa = True, progressbar=False):
    '''
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
    INPUT
        data: image stack, must be 2D, must be in real space
        mode: whether to apply alignment parameters. Default mode='a' means apply parameters
        rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. This is only relevant for the case where images are 2D, in which case rot_method can be either rot_shift2D or rotshift2dg, with the default being rot_shift2D. If images are 3D, rot_shift3D will be used to rotate/shift the images.
        interp: interpolation method to use for rot_method when applying alignment parameters.
        i1: index of first image to be used.
        i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
        use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
        use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
        snr: signal to noise ratio, default 1.0
    
    OUTPUT
        tavg: The best estimate (Wiener filter) given the image series and estimated CTF parms, in real space.
        var: Variance (in real space) calculated as follows, [1/(n-1)]*[sum_j { F[(H_j*(O_j - F^{-1}(H_j*tavg))^2]/SUM_CTF^2} where O_j is the j-th image in real space, F^{-1} denotes inverse fourier transform operator, and H_j is the CTF of the j-th image
    
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
        ERROR("images must be 2D for CTF correction.....exiting","avgvar_ctf",1)

    if img.get_attr_default('ctf_applied', 0) == 1:
        ERROR("data cannot be ctf-applied....exiting","avgvar_ctf",1)

    if inmem:
        data_nima = len(data)
    else:
        data_nima = EMUtil.get_image_count(data)

    if i2 == 0: i2 = data_nima-1
    if dopa:
        nx2 = nx*2
        ny2 = ny*2
    else:
        nx2 = nx
        ny2 = ny
    ave = EMAN2.EMData(nx2, ny2, 1, False)
    ctf_2_sum = EMAN2.EMData(nx2, ny2, 1, False)
    nima = 0
    for i in tqdm.tqdm(range(i1, i2+1), unit='img', disable=(not progressbar), file=sys.stdout):
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

        img = sp_utilities.pad(img, nx2, ny2, 1, background = "circumference")
        sp_fundamentals.fftip(img)
        EMAN2_cppwrap.Util.add_img(ave, sp_filter.filt_ctf(img, ctf_params))
        EMAN2_cppwrap.Util.add_img2(ctf_2_sum, sp_morphology.ctf_img(nx2, ctf_params))

    ctf_2_sum += 1.0/snr
    EMAN2_cppwrap.Util.div_filter(ave, ctf_2_sum)

    # calculate variance in real space
    tvar = sp_utilities.model_blank(nx, ny, nz)
    for i in tqdm.tqdm(range(i1, i2+1), unit='img', disable=(not progressbar), file=sys.stdout):
        if not(use_odd) and i%2 == 1: continue
        if not(use_even) and i%2 == 0: continue
        if inmem: img = data[i].copy()
        else: img = sp_utilities.get_im(data, i)

        ctf_params = img.get_attr("ctf")

        if (mode == 'a'):
            angle, sx, sy, mirror, scale = sp_utilities.get_params2D(img)
            img = sp_fundamentals.rot_shift2D(img, angle, sx, sy, mirror, scale, interp)
            ctf_params.dfang += angle
            if mirror == 1:  ctf_params.dfang = 270.0-ctf_params.dfang

        img = sp_utilities.pad(img, nx2, ny2, 1, background = "circumference")
        sp_fundamentals.fftip(img)
        img = img-sp_filter.filt_ctf(ave, ctf_params, dopa)
        #Util.div_filter(img, ctf_2_sum)
        img = sp_fundamentals.window2d(sp_fundamentals.fft(img),nx,ny)
        #Util.add_img(totv, img)
        EMAN2_cppwrap.Util.add_img2(tvar, img)
    #Util.mul_scalar(tvar, float(nima)*nima/(nima-1)) # the strange factor is due to the fact that division by ctf^2 is equivalent to division by nima
    EMAN2_cppwrap.Util.mul_scalar(tvar, 1.0/float(nima))
    return  sp_fundamentals.window2d(sp_fundamentals.fft(ave),nx,ny) , tvar#,(tvar - totv*totv/nima), tvar, totv,tavg

def avgvar_flipctf(data, mode='a', interp='quadratic', i1=0, i2=0, use_odd=True, use_even=True, snr=1.0, dopa = True, progressbar=False):
    '''
    Averages data while applying phase-flip.
    Adapted from sp_statistics.avgvar_ctf
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
    Parameters
        data: image stack, must be 2D, must be in real space
        mode: whether to apply alignment parameters. Default mode='a' means apply parameters
        rot_method: specifies the function by which images are rotated/shifted if alignment parameters are to be applied. This is only relevant for the case where images are 2D, in which case rot_method can be either rot_shift2D or combined2dg, with the default being rot_shift2D. If images are 3D, rot_shift3D will be used to rotate/shift the images.
        interp: interpolation method to use for rot_method when applying alignment parameters.
        i1: index of first image to be used.
        i2: index of last image to be used. If i2 = 0, then i2 defaults to one less than number of images in the data
        use_odd: images with indices between i1 and i2 which are odd are used if and only if use_odd is set to True. Default is True.
        use_even: images with indices between i1 and i2 which are even are used if and only if use_even is set to True. Default is True.
        snr: signal to noise ratio, default 1.0
        dopa: flag to pad images by 2
        progressbar: flag to show progress bar
    
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
        data_nima = EMUtil.get_image_count(data)

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
    for i in tqdm.tqdm(range(i1, i2+1), unit='img', disable=(not progressbar), file=sys.stdout):
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

        img = sp_utilities.pad(img, nx2, ny2, 1, background = "circumference")
        sp_fundamentals.fftip(img)
        img= sp_filter.filt_ctf(img, ctf_params, binary=1)
        
        # Update average and variance
        EMAN2_cppwrap.Util.add_img(ave, img)
        EMAN2_cppwrap.Util.add_img2( var, sp_fundamentals.fft(img) )

    EMAN2_cppwrap.Util.mul_scalar(ave, 1.0 /float(nima) )
    avg_padded= sp_fundamentals.fft(ave)
    var_padded= (var - avg_padded*avg_padded*nima)/(nima-1)
    
    # Back to original size
    avg_unpad= sp_fundamentals.window2d(avg_padded, nx, ny)
    var_unpad= sp_fundamentals.window2d(var_padded, nx, ny)
    
    return avg_unpad, var_unpad, ctf_params.apix
    
def punish_stack(input_stack, output_stack=None, outbdb=None, idim0=-1, first=0, last=-1, 
                ctf_profile_doc=None, plotfile=None, log=None, verbosity=2, progressbar=True, 
                debug=False, padfile=None, straight_ctf_doc=None, straight_ctf_img=None, 
                flip_ctf_doc=None, flip_ctf_img=None, trap_ctf_doc=None, trap_ctf_img=None):
    """
    Multiples image stack by CTF but while maintaining amplitudes until the first extremum.
    Adapted from sp_ctf_punisher.py
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
    Arguments:
        input_stack : Input stack
        output_stack : Output stack filename, e.g., MRCS, HDF, etc.
        outbdb : (optional) BDB to which header information will be written
        idim0 : Unpadded image dimension
        first : First image in stack to process
        last : Last image in stack to process
        ctf_profile_doc : Doc with CTF profiles
        plotfile : CTF profile plot file
        log : Logger object
        verbosity : How much additional information to write to screen (0..7)
        progressbar : (boolean) Flag to show progress bar
        debug : (boolean) Write example of images: unfiltered, CTF-multiplied, phase-flipped, and trapped CTF
        padfile : Padded image example (debug mode)
        straight_ctf_doc : CTF profile doc (debug mode)
        straight_ctf_img : CTF-multiplied image example (debug mode)
        flip_ctf_doc : Phase-flipping 1D profile doc (debug mode)
        flip_ctf_img : Phase-flipped image example (debug mode)
        trap_ctf_doc : Trapped-CTF 1D profile doc (debug mode)
        trap_ctf_img : Trapped-CTF image example (debug mode)
        
    Returns:
        List of EMData objects
    """
    
    # Initialize
    partnum= 0
    
    # Check whether BDB or list
    if type(input_stack) == type(""):
        inmem = False
    else:
        inmem = True
    
    # Get dimensions
    if idim0<=0: 
        if inmem:
            idim0= img = input_stack[0]['nx']
        else:
            idim0= sp_utilities.get_im(input_stack, 0)['nx']
    
    # Get last image if not specified
    if last<=0: 
        if inmem:
            last= len(input_stack)
        else:
            last= EMUtil.get_image_count(input_stack)
    
    assert first<last, "ERROR!! Image range %s < last" % (first, last)
    
    disableTF= verbosity<1 or verbosity>=3
    if not progressbar : disableTF= True
    
    # If we pretend the BDB file is in a subdirectory by substituting a '/' for the "#", then makerelpath is easy
    if outbdb: 
        assert outbdb[:4]=='bdb:', "ERROR!! BDB '%s' does not start with 'bdb:'" % outbdb
        output_bdb_stem= outbdb[4:].replace("#",'/')
        new_bdb_dict= db_open_dict(outbdb)
    
    image_list= []
    
    # Loop through images
    for idx in tqdm.tqdm(range(first, last), unit='img', disable=disableTF, file=sys.stdout):
        # Read image
        if inmem:
            img_obj= input_stack[idx]
        else:
            img_obj= sp_utilities.get_im(input_stack, idx)
        
        # If first image, get image dimension
        if idx==first:
            idim0= img_obj['nx']
            
            # Will pad image by 2
            paddim= 2*idim0
        
        # Pad image by 2
        padded_img= padimgby2(img_obj)
        
        # Read CTF parameters (If sign=1 (default), starts out positive.)
        ctf_params= img_obj['ctf']
        if verbosity >= 4: 
            print_log_msg("Particle %s: Defocus %s" % (idx, ctf_params.to_dict()['defocus']) , log)
            
        # Generate 1D CTF profile
        straight_ctf= sp_morphology.ctf_1d(paddim, ctf_params, sign=1)
        
        # List has sqrt(2)*nx elements
        straight_firstn= straight_ctf[:idim0]
        
        # Apply straight CTF
        ctfd_img= sp_filter.filt_table(padded_img, straight_ctf)
        
        # Set the CTF to 1 until the first maximum
        punished_ctf, flipped_ctf, combined_ctfs, _ = trapctf(copy.copy(straight_ctf), log=log, verbosity=verbosity)
        # (Weird things happen if I don't send a copy of the list.)
        
        # Apply punished CTF
        punished_img= sp_filter.filt_table(padded_img, punished_ctf)
            
        # Write profile to disk
        if ctf_profile_doc != None and idx == 0 : sp_utilities.write_text_row(combined_ctfs, ctf_profile_doc)
        
        # Write diagnostic images
        if debug and idx == 0 : 
            padded_img.write_image(padfile)
            sp_utilities.write_text_row(straight_ctf, straight_ctf_doc)
            ctfd_img.write_image(straight_ctf_img)
            sp_utilities.write_text_row(flipped_ctf, flip_ctf_doc)
            
            #flipped_img= filt_ctf(padded_img, ctf_params, dopad=False, binary=1)  
            ## (This looks a little different, but using filt_table for consistency.)
            
            flipped_img= sp_filter.filt_table(padded_img, flipped_ctf)
            flipped_img.write_image(flip_ctf_img)
            
            sp_utilities.write_text_row(punished_ctf, trap_ctf_doc)
            punished_img.write_image(trap_ctf_img)
            
        # Crop back down to original size
        cropped_img= sp_fundamentals.window2d(punished_img, idim0, idim0)
        
        # Get last image number
        if output_stack:
            if os.path.exists(output_stack):
                next_imgnum= EMUtil.get_image_count(output_stack)
            else:
                next_imgnum= 0
            cropped_img.write_image(output_stack, next_imgnum)
        
        image_list.append(cropped_img)
        
        # If output BDB supplied, write to BDB
        if outbdb:
            img_dict= img_obj.get_attr_dict()
            assert output_stack != None, "ERROR!! Output stack name is required if BDB is written"
            img_dict["data_path"]= makerelpath(output_bdb_stem, output_stack)
            img_dict["ptcl_source_coord_id"]= idx - first
            new_bdb_dict[idx]= img_dict
        
        # Plot
        if idx == 0 and plotfile != None:
            punished_firstn= punished_ctf[:idim0]
            flipped_firstn= flipped_ctf[:idim0]
            
            try:
                plt.plot(range(idim0), straight_firstn, color='lime', label='Straight CTF')
                plt.plot(range(idim0), flipped_firstn, color='blue', label='Flipped CTF')
                plt.plot(range(idim0), punished_firstn, color='orchid', label='Punished CTF')
                plt.legend(loc="upper right")
                plt.ylim(-1.1, 1.6)
                plt.xlabel('Fourier radius, pixels')
                plt.savefig(plotfile)
            
            except RuntimeError:
                print_log_msg("WARNING! Can't create CTF plot (a machine-specific problem). Skipping...", log)
    # End image loop
    
    if outbdb : db_close_dict(outbdb)
    
    if output_stack and verbosity>=3:
        mic_base= os.path.basename(output_stack)
        print_log_msg('%s: first %s, last %s' % (mic_base, first, last-1), log)
        
    return image_list

def padimgby2(input_img_obj, factor=2, verbosity=0, log=None):
    """
    Pads image. 
    Original image will be centered in the padded one.
    
    Copied from sp_center2d3d.py
    Adapted from sp_ctf_punisher.py
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
    Arguments:
        input_img_obj : input EMData object
        factor : (int) multiplication factor
    Returns:
        padimg : padded image
    """
    
    # Assume square
    idim = input_img_obj['nx']
    if verbosity >= 7 : print_log_msg('Starting dimension: %s' % idim, log)

    # Pad by specified padding factor
    paddim = idim*factor
    if verbosity >= 7 : print_log_msg('Computed padded dimension: %s' % paddim, log)

    # Pad (will be centered in new box)
    padimg = sp_utilities.pad(input_img_obj, paddim, paddim, background='circumference')
    if verbosity >= 7 : print_log_msg('Confirmed dimension: %s' % padimg['nx'], log)
    
    return padimg
    
def trapctf(ctflist, log=None, verbosity=0):
    """
    Sets the CTF to 1 until the first maximum.
    
    Adapted from SPIDER's trapctf.cor
    Adapted from sp_ctf_punisher.py
    Modified 2020-07-03 for SPHIRE 1.4 (Python 3)
    
    Arguments:
        ctflist : List of CTF values, assumed to start out positive, i.e. ctf_1d(sign=1)
    """
    
    # Make sure the CTF increases
    assert ctflist[1] > 0, "Uh oh!! CTF (%s) is negative!" % ctflist[1]
    
    num_radii= len(ctflist)
    
    # initialize
    firstmax= num_radii + 1  
    max_ctf= -2
    foundMax= False
    trapctf_list= []
    flipctf_list= []
    combined_list= []
    
    # Loop though list to look for maximum
    for idx in range(num_radii):
        curr_straightctf= ctflist[idx]
        
        # Phase-flipped profile
        if curr_straightctf != 0.0:
            curr_flippedctf= abs(curr_straightctf)/curr_straightctf
        else:
            curr_flippedctf= 0.0
        flipctf_list.append(curr_flippedctf)
        
        if foundMax == False and curr_straightctf > max_ctf:
            max_ctf= curr_straightctf
            
            # Overwrite ctflist
            ctflist[idx]= 1.0
            
            # Will also return after & before list
            curr_trappedctf= 1.0
        else:
            curr_trappedctf= curr_straightctf
            
            if foundMax == False:
                foundMax= True
                firstmaxradius= idx - 1
        trapctf_list.append(curr_trappedctf)
        
        combined_list.append([curr_trappedctf, curr_straightctf, curr_flippedctf])
        
        if verbosity >= 6 : 
            print_log_msg("radius %s: %s, foundMax %s" % (idx, combined_list[-1], foundMax), log)
    # End first-maximum loop
    
    if verbosity >= 5 : print_log_msg('First max at %s Fpx' % firstmaxradius, log)
    
    return trapctf_list, flipctf_list, combined_list, firstmaxradius
    
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
        'input_star_file', 
        type=str, 
        help='Input STAR file')
    
    parser.add_argument(
        'outdir',
        type=str, 
        help='Output directory')
    
    parser.add_argument(
        '--stackdir', 
        type=str, 
        help='Directory where particle stacks are located relative to path in STAR file')
    
    parser.add_argument(
        '--ctf', 
        type=str, 
        default=None, 
        help='CTF-correction options for averages: flip, wiener, punish')
    
    parser.add_argument(
        '--unaligned', 
        type=str, 
        default=None, 
        help='STAR file for unaligned images (for debugging)')
    
    parser.add_argument('--negate', 
        action="store_true", 
        help='Negate signs for shifts from STAR file')
    
    parser.add_argument('--verbosity', "-v", 
        type=int, 
        default=3, 
        help='Increase verbosity (0..5)')
    
    return parser.parse_args()

def main():
    options = parse_command_line()
    
    ##print args, options  # (Everything is in options.)
    #print(options)
    ##print('LOGFILE',sp_global_def.LOGFILE)
    #exit()

    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
    if RUNNING_UNDER_MPI:
            mpi.mpi_init( 0, [] )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp( "Start" )	

    relion_avg(options.input_star_file, outdir=options.outdir, options=options, verbosity=options.verbosity, usempi=RUNNING_UNDER_MPI)

    sp_global_def.print_timestamp( "Finish" )
    if RUNNING_UNDER_MPI:
            mpi.mpi_finalize()

if __name__ == "__main__":
    main()
