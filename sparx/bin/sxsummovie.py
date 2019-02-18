#! /usr/bin/env python
from __future__ import print_function
#
# Author 2016  Markus Stabrin (markus.stabrin@mpi-dortmund.mpg.de)
# Copyright (C) 2019 Max planck institute for molecular physiology, Dortmund
#

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

from sys import argv
from os import mkdir, path, system, remove, rmdir
from glob import glob
from numpy import genfromtxt
import numpy
import time
import subprocess
import global_def
from global_def import SPARXVERSION, ERROR
from optparse import OptionParser, SUPPRESS_HELP
from utilities import create_summovie_command


def main():

    # Parse the Options
    progname = path.basename(argv[0])
    usage = progname + """ summovie_path input_micrograph_pattern input_shift_pattern output_directory
    --selection_list
    --nr_frames=nr_frames
    --first
    --last
    --pixel_size=pixel_size
    --nr_threads
    --apply_dose_filter
    --voltage=voltage
    --exposure_per_frame=exposure_per_frame
    --pre_exposure=pre_exposure
    --dont_restore_noise

    sxsummovie exists only in non-MPI version.

    Perform summovie without dose filtering.

    sxsummovie.py ~/my_app/summovie 'outdir_unblur/corrsum/micrograph_*_frames_sum.mrc' 'outdir_unblur/shift/micrograph_*_frames_shift.txt'
    outdir_summovie --nr_frames=24 --pixel_size=1.19 --nr_threads=1

    Perform summovie without dose filtering and with less frames.

    sxsummovie.py ~/my_app/summovie 'outdir_unblur/corrsum/micrograph_*_frames_sum.mrc' 'outdir_unblur/shift/micrograph_*_frames_shift.txt'
    outdir_summovie --nr_frames=24 --first=3 --last=15 --pixel_size=1.19 --nr_threads=1

    Perform summovie with dose filtering and with less frames.

    sxsummovie.py ~/my_app/summovie 'outdir_unblur/corrsum/micrograph_*_frames_sum.mrc' 'outdir_unblur/shift/micrograph_*_frames_shift.txt'
    outdir_summovie --nr_frames=24 --first=3 --last=15 --pixel_size=1.19 --nr_threads=1 --apply_dose_filter --voltage=300 --exposure_per_frame=2 --pre_exposure=0
    """

    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option('--selection_list',     type='str',          default='',        help='Micrograph selecting list (SPHIRE specific): Specify a name of micrograph selection list text file. The file extension must be \'.txt\'. If this is not provided, all files matched with the micrograph name pattern will be processed. (default none)')
    parser.add_option('--nr_frames',          type='int',          default=3,         help='Number of movie frames: The number of movie frames in each input micrograph. (default 3)')
    parser.add_option('--first',              type='int',          default=1,         help='First movie frame: First movie frame for summing. (default 1)')
    parser.add_option('--last',               type='int',          default=-1,        help='Last movie frame: Last movie frame for summing. (default -1)')
    parser.add_option('--sum_suffix',         type='str',          default='_sum',    help=SUPPRESS_HELP)
    parser.add_option('--pixel_size',         type='float',        default=-1.0,      help='Pixel size [A]: The pixel size of input micrographs. (default required float)')
    parser.add_option('--nr_threads',         type='int',          default=1,         help='Number of threads: The number of threads summovie can use. The higher the faster, but it requires larger memory. (default 1)')
    parser.add_option('--apply_dose_filter',  action='store_true', default=False,     help='Apply dose filter step: Requires voltage, exposure per frame, and pre exposure options. (default False)')
    parser.add_option('--voltage',            type='float',        default=300.0,     help='Microscope voltage (dose filter) [kV]: The acceleration voltage of microscope used for imaging. (default 300.0)')
    parser.add_option('--exposure_per_frame', type='float',        default=2.0,       help='Per frame exposure (dose filter) [e/A^2]: The electron dose per frame in e/A^2. (default 2.0)')
    parser.add_option('--pre_exposure',       type='float',        default=0.0,       help='Pre-exposure (dose filter) [e/A^2]: The electron does in e/A^2 used for exposure prior to imaging .(default 0.0)')
    parser.add_option('--frc_suffix',         type='string',       default='_frc',    help=SUPPRESS_HELP)
    parser.add_option('--dont_restore_noise', action='store_true', default=False,     help='Do not restore noise power: Do not restore noise power. (default False)')
    parser.add_option('--summovie_ready',     action='store_true', default=False,     help=SUPPRESS_HELP)

    # list of the options and the arguments
    (options, args) = parser.parse_args(argv[1:])

    global_def.BATCH = True

    # If there arent enough arguments, stop the script
    if len(args) != 4:
        ERROR("see usage " + usage, 1)

    # Convert the realtive parts to absolute ones
    summovie_path = path.realpath(args[0]) # summovie_path
    input_image = path.realpath(args[1])   # input_micrograph_pattern
    input_shift = path.realpath(args[2])   # input_shift_pattern
    output_dir = path.realpath(args[3])    # output_directory

    # If the summovie executable file does not exists, stop the script
    if not path.exists(summovie_path):
        ERROR(
            'Summovie directory does not exist, please change' +
            ' the name and restart the program.', 'sxsummovie.py', 1
            )

    # If the output directory exists, stop the script
    if path.exists(output_dir):
        ERROR(
            'Output directory exists, please change' +
            ' the name and restart the program.', 'sxsummovie.py', 1
            )

    # If the input file does not exists, stop the script
    file_list = glob(input_image)
    shift_list = glob(input_shift)

    if not file_list:
        ERROR(
            'Input micrograph file(s) does not exist, please change' +
            ' the name and restart the program.', 'sxsummovie.py', 1
            )

    if not shift_list:
        ERROR(
            'Input shift file(s) does not exist, please change' +
            ' the name and restart the program.', 'sxsummovie.py', 1
            )

    # Output paths
    if options.apply_dose_filter:
        output_path = '{:s}/corrsum_dose_filtered'.format(output_dir)
    else:
        output_path = '{:s}/corrsum'.format(output_dir)
    frc_path = '{:s}/frc'.format(output_dir)
    log_path = '{:s}/logfiles'.format(output_dir)
    temp_path = '{0}/temp'.format(output_dir)

    # Split the path of the image name at the "/" Characters.
    # The last entry contains the micrograph name.
    # Split the micrograph name at the wildcard character for the
    # prefix and suffix.
    input_mic_split = input_image.split('/')
    input_mic_name = input_mic_split[-1].split('*')
    input_shift_name = input_shift.split('*')

    if len(input_mic_name) != 2 or len(input_shift_name) != 2:
        ERROR(
            'Too many wildcard arguments.' +
            'Please use exactly one * in the pattern.', 'sxsummovie.py',
            1
            )

    # Get the input directory
    if len(input_mic_split) != 1:
        input_dir = input_image[:-len(input_mic_split[-1])]
    else:
        input_dir = ''

    # Create output directorys
    if not path.exists(output_dir):
        mkdir(output_dir)
    if not path.exists(output_path):
        mkdir(output_path)
    if not path.exists(frc_path):
        mkdir(frc_path)
    if not path.exists(temp_path) and not options.summovie_ready:
        mkdir(temp_path)
    if not path.exists(log_path):
        mkdir(log_path)

    # shift wildcard list
    shift_wildcard = [entry[len(input_shift_name[0]):-len(input_shift_name[-1])] \
            for entry in shift_list]

    # Just use shifts that have a micrograph and vise versa
    mic_list = [
            entry for entry in file_list \
                if entry[len(input_dir) + len(input_mic_name[0]):-len(input_mic_name[-1])]\
                in shift_wildcard] 

    # If micrograph list is provided just process the images in the list
    selection_file = options.selection_list
    if selection_file:
        # Import list file
        try:
            selection = genfromtxt(selection_file, dtype=None)
        except TypeError:
            ERROR('no entrys in micrograph list file {0}'.format(selection_file), 'sxsummovie.py', 1)
        # List of files which are in pattern and list
        mic_list = [
                entry for entry in mic_list \
                if entry[len(input_dir):] in selection and \
                path.exists(entry)
                ]
        # If no match is there abort
        if len(mic_list) == 0:
            ERROR(
                'no files in {0} matched the micrograph file pattern:\n'.format(selection_file), 'sxsummovie.py',
                1
                )

    option_dict = {
        'summovie_path': summovie_path,
        'mic_list': mic_list,
        'mic_prefix': input_mic_name[0],
        'mic_suffix': input_mic_name[1],
        'shift_prefix': input_shift_name[0],
        'shift_suffix': input_shift_name[1],
        'input_dir': input_dir,
        'output_dir': output_dir,
        'output_path': output_path,
        'frc_path': frc_path,
        'log_path': log_path,
        'temp_path': temp_path,
        'nr_frames': options.nr_frames,
        'sum_suffix': options.sum_suffix,
        'first': options.first,
        'last':options.last,
        'pixel_size': options.pixel_size,
        'apply_dose_filter': options.apply_dose_filter,
        'exposure_per_frame': options.exposure_per_frame,
        'voltage': options.voltage,
        'pre_exposure': options.pre_exposure,
        'frc_suffix': options.frc_suffix,
        'dont_restore_noise': options.dont_restore_noise,
        'nr_threads': options.nr_threads,
        'summovie_ready': options.summovie_ready,
        'verbose': True
        }

    # Run summovie
    run_summovie(
        opt=option_dict
        )

    if not options.summovie_ready:
        # Remove temp folder
        for entry in glob('{0}/*'.format(temp_path)):
            remove(entry)
        rmdir(temp_path)

    print('All Done!')

    global_def.BATCH = False


def run_summovie(
        opt
        ):

    # Lists to write the text files later
    micrograph_list = []

    # Get the number of files
    nr_files = len(opt['mic_list'])

    # Timeing stuff
    time_start = time.time()
    time_list = []

    # Loop over all files
    for index, inputfile in enumerate(sorted(opt['mic_list'])):

        # Check, if there is an prefix and suffix.
        # If there is more then one entry: the suffix is the last one.
        # Otherwhise its just the one after the dot.
        input_suffix = inputfile.split('/')[-1].split('.')[-1]
        # First output to introduce the programm
        if opt['verbose'] and index == 0:
            print(
                    'Progress: 0.0%;  Time: --h:--m:--s/--h:--m:--s;  Summovie started!'
                )

        # Time begin
        t1 = time.time()

        # Get the output names
        file_name = inputfile[len(opt['input_dir']):-len(opt['mic_suffix'])]
        file_wildcard = file_name[
                len(opt['mic_prefix']):
                ]
        micrograph_name = '{0}/{1}{2}.mrc'.format(
            opt['output_path'], file_name, opt['sum_suffix']
            )
        frc_name = '{0}/{1}{2}.txt'.format(
                opt['frc_path'], file_name, opt['frc_suffix']
                )
        shift_name = '{0}{1}{2}'.format(
                opt['shift_prefix'], file_wildcard, opt['shift_suffix']
                )
        if not opt['summovie_ready']:
            temp_name = '{0}/{1}{2}.mrc'.format(
                    opt['temp_path'], file_name, opt['sum_suffix']
                    )
        else:
            temp_name = inputfile
        log_name = '{0}/{1}.log'.format(
                opt['log_path'], file_name
                )
        error_name = '{0}/{1}.err'.format(
                opt['log_path'], file_name
                )
        # Append the names to the lists
        micrograph_list.append('{0}{1}.mrc'.format(
            file_name, opt['sum_suffix'])
            )

        # First build the summovie command
        summovie_command = create_summovie_command(
            temp_name,
            micrograph_name,
            shift_name,
            frc_name,
            opt
            )

        # Export the number of threads
        export_threads_command = []

        # Export
        export_threads_command.append('export')
        # Nr of threads
        export_threads_command.append('OMP_NUM_THREADS={0}'.format(
            opt['nr_threads']
            ))

        if not opt['summovie_ready']:
            # Do a e2proc3d.py
            e2proc3d_command = []

            # e2proc3d
            e2proc3d_command.append('e2proc3d.py')
            # inputfile
            e2proc3d_command.append('{0}'.format(inputfile))
            # outputfile
            e2proc3d_command.append('{0}'.format(temp_name))


        # Translate the command to single strings
        if not opt['summovie_ready']:
            e2proc3d_command = r' '.join(e2proc3d_command)
        export_threads_command = r' '.join(export_threads_command)
        summovie_command = '\n'.join(summovie_command)

        # Build full command
        if not opt['summovie_ready']:
            full_command = r'{0}; {1}; echo "{2}" | {3}'.format(
                export_threads_command,
                e2proc3d_command,
                summovie_command,
                opt['summovie_path']
                )
        else:
            full_command = r'{0}; echo "{1}" | {2}'.format(
                export_threads_command,
                summovie_command,
                opt['summovie_path']
                )

        # Remove temp summovie files
        temp_summovie_files = glob('.SumMovie*')
        for entry in temp_summovie_files:
            remove(entry)

        with open(log_name, 'w') as f:
            with open(error_name, 'w') as e:
                # Execute Command
                subprocess.Popen(
                    [full_command], shell=True,
                    stdout=f,
                    stderr=e
                    ).wait()

        # Remove temp summovie files
        temp_summovie_files = glob('.SumMovie*')
        for entry in temp_summovie_files:
            remove(entry)
        if not opt['summovie_ready']:
            if path.exists(temp_name):
                # Remove temp file
                remove(temp_name)
            else:
                ERROR('e2proc2d.py error. File was not created:\n{0}'.format(inputfile), 'sxsummovie.py', 0)

        time_list.append(time.time() - t1)
        
        # Check if SumMovie finished cleanly
        with open(log_name, 'r') as r:
            clean = False
            for line in r:
                if 'SumMovie finished cleanly.' in line:
                    clean = True
                    break
        if clean:
            print('SumMovie finished cleanly.')
        else:
            ERROR(
                'sum movie error. check the logfile for more information: {0}'.format(
                    log_name), 'sxsummovie.py', 0
                )

        # Do progress output
        if opt['verbose']:
            percent = round(100 * (index + 1) / float(nr_files), 2)
            estimated_time = \
                nr_files * sum(time_list) / float(len(time_list))
            estimated_time_h = estimated_time // 3600
            estimated_time_m = (estimated_time - estimated_time_h*3600) // 60
            estimated_time_s = (
                    estimated_time -
                    estimated_time_h*3600 -
                    estimated_time_m*60
                    )
            current_time = time.time() - time_start
            current_time_h = current_time // 3600
            current_time_m = (current_time - current_time_h*3600) // 60
            current_time_s = (
                    current_time -
                    current_time_h*3600 -
                    current_time_m*60
                    )
            print((
                'Progress: {0:.2f}%;  Time: {1:.0f}h:{2:.0f}m:{3:.0f}s/{4:.0f}h:{5:.0f}m:{6:.0f}s;  Micrograph done:{7}'.format(
                    percent,
                    current_time_h,
                    current_time_m,
                    current_time_s,
                    estimated_time_h,
                    estimated_time_m,
                    estimated_time_s,
                    file_name
                    )
                ))


    # Write micrograph list
    with open('{0}/summovie_micrographs.txt'.format(opt['output_dir']), 'w') as f:
        for entry in sorted(micrograph_list):
            f.write('{0}\n'.format(entry))


if __name__ == '__main__':
    main()
