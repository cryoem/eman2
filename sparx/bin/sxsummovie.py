#! /usr/bin/env python
#
# Copyright (C) 2016  Markus Stabrin (markus.stabrin@mpi-dortmund.mpg.de)
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


def main():

    # Parse the Options
    progname = path.basename(argv[0])
    usage = progname + """ summovie input_image input_shift output
    --input_micrograph_list_file
    --input_shift_list_file
    --nr_frames=nr_frames
    --pixel_size=pixel_size
    --first
    --last
    --apply_dose_filter
    --exposure_per_frame=exposure_per_frame
    --voltage=voltage
    --pre_exposure=pre_exposure
    --dont_restore_noise
    --nr_threads'

    sxsummovie exists in non-MPI version.

    """

    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option('--selection_list',         type='str',          default='',    help='not a summovie option: input selection micrograph list file: Extension of input micrograph list file must be ".txt". If this is not provided, all files matched with the micrograph name pattern will be processed. (default none)')
    parser.add_option('--nr_frames',          type='int',          default=3,         help='number of frames in the set of micrographs')
    parser.add_option('--sum_suffix',         type='str',          default='_sum',    help=SUPPRESS_HELP)
    parser.add_option('--pixel_size',         type='float',        default=-1.0,      help='pixel size [A]')
    parser.add_option('--apply_dose_filter',        action='store_true', default=False,     help='apply dose filter options')
    parser.add_option('--exposure_per_frame', type='float',        default=2.0,       help='exposure per frame [e/A^2]')
    parser.add_option('--voltage',            type='float',        default=300.0,     help='accelerate voltage [kV]')
    parser.add_option('--pre_exposure',       type='float',        default=0.0,       help='pre exposure amount [e/A^2]')
    parser.add_option('--frc_suffix',         type='string',       default='_frc',    help=SUPPRESS_HELP)
    parser.add_option('--dont_restore_noise',      action='store_true', default=False,     help='do not restore noise power')
    parser.add_option('--nr_threads',         type='int',          default=1,         help='number of threads summovie is allowed to use. The higher the faster, but it also needs a higher amount of memory.')
    parser.add_option('--summovie_ready',        action='store_true', default=False,      help=SUPPRESS_HELP)

    # list of the options and the arguments
    (options, args) = parser.parse_args(argv[1:])

    global_def.BATCH = True

    # If there arent enough arguments, stop the script
    if len(args) != 4:
        ERROR("see usage " + usage, 1)

    # Convert the realtive parts to absolute ones
    summovie_path = path.realpath(args[0])
    input_image = path.realpath(args[1])
    input_shift = path.realpath(args[2])
    output_dir = path.realpath(args[3])

    # If the summovie executable file does not exists, stop the script
    if not path.exists(summovie_path):
        ERROR(
            'Summovie directory does not exist, please change' +
            ' the name and restart the program.', 1
            )

    # If the output directory exists, stop the script
    if path.exists(output_dir):
        ERROR(
            'Output directory exists, please change' +
            ' the name and restart the program.', 1
            )

    # If the input file does not exists, stop the script
    file_list = glob(input_image)
    shift_list = glob(input_shift)

    if not file_list:
        ERROR(
            'Input micrograph file(s) does not exist, please change' +
            ' the name and restart the program.', 1
            )

    if not shift_list:
        ERROR(
            'Input shift file(s) does not exist, please change' +
            ' the name and restart the program.', 1
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

    # Get the input directory
    if len(input_split) != 1:
        input_dir = input_image[:-len(input_split[-1])]
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
    shift_wildcard = [entry for entry in 

    # Just use shifts that have a micrograph and vise versa
    mic_list = [
            entry for entry in mic_list \
                if entry[len(input_dir) + len(input_mic_name[0]):-len(input_muc_name[-1])]\
                in 

    # If micrograph list is provided just process the images in the list
    mic_file = options.selection_list
    if mic_file:
        # Import list file
        try:
            selection = genfromtxt(mic_file, dtype=None)
        except TypeError:
            ERROR('no entrys in micrograph list file {0}'.format(mic_file), 1)
        # List of files which are in pattern and list
        mic_file_list = [
                entry for entry in file_list \
                if entry[len(input_dir):] in selection and \
                path.exists(entry)
                ]
        # If no match is there abort
        if len(file_list) == 0:
            ERROR(
                'no files in {0} matched the micrograph file pattern:\n'.format(mic_file),
                1
                )

    option_dict = {
        'summovie_path': summovie_path,
        'input_image': input_image,
        'input_dir': input_dir,
        'input_shift': input_shift,
        'output_dir': output_dir,
        'output_path': output_path,
        'frc_path': frc_path,
        'temp_path': temp_path,
        'log_path': log_path,
        'file_list': file_list,
        'first': options.first,
        'last': options.last,
        'pixel_size': options.pixel_size,
        'nr_frames': options.nr_frames,
        'sum_suffix': options.sum_suffix,
        'apply_dose_filter': options.apply_dose_filter,
        'voltage': options.voltage,
        'pre_exposure': options.pre_exposure,
        'frc_suffix': options.frc_suffix,
        'dont_restore_noise': options.dont_restore_noise,
        'nr_threads': options.nr_threads,
        'summovie_ready': options.summovie_ready

    # Run summovie
    run_summovie(
        options=option_dict
        )

    if not options.summovie_ready:
        # Remove temp folder
        for entry in glob('{0}/*'.format(temp_path)):
            remove(entry)
        rmdir(temp_path)

    print('All Done!')

    global_def.BATCH = False


def prepare_summovie(
        summovie_path,
        input_image,
        input_dir,
        output_dir,
        corrected_path,
        uncorrected_path,
        shift_path,
        frc_path,
        temp_path,
        log_path,
        file_list,
        options
        ):

    # Lists to write the text files later
    micrograph_list = []
    shift_list = []
    if options.save_frames:
        frames_list = []

    # If micrograph list is provided just process the images in the list
    mic_list = options.input_micrograph_list_file
    if mic_list:
        # Import list file
        try:
            set_selection = genfromtxt(mic_list, dtype=None)
        except TypeError:
            ERROR('no entrys in list file {0}'.format(mic_list), 1)
        # List of files which are in pattern and list
        file_list = [
                entry for entry in file_list \
                if entry[len(input_dir):] in set_selection and \
                path.exists(entry)
                ]
        # If no match is there abort
        if len(file_list) == 0:
            ERROR(
                'no files in {0} matched the file pattern:\n'.format(mic_list),
                1
                )
    # Get the number of files
    nr_files = len(file_list)

    # Timeing stuff
    time_start = time.time()
    time_list = []

    # Loop over all files
    for index, inputfile in enumerate(sorted(file_list)):

        # Check, if there is an prefix and suffix.
        # If there is more then one entry: the suffix is the last one.
        # Otherwhise its just the one after the dot.
        input_suffix = inputfile.split('/')[-1].split('.')[-1]
        # First output to introduce the programm
        if index == 0:
            print(
                    'Progress: 0.0%;  Time: --h:--m:--s/--h:--m:--s;  Unblur started!'
                )

        # Time begin
        t1 = time.time()

        # Get the output names
        file_name = inputfile[len(input_dir):-len(input_suffix) - 1]
        if options.skip_dose_filter:
            micrograph_name = '{0}/{1}{2}.mrc'.format(
                    uncorrected_path, file_name, options.sum_suffix
                    )
            frames_name = '{0}/{1}{2}.mrc'.format(
                    uncorrected_path, file_name, options.frames_suffix
                    )
        else:
            micrograph_name = '{0}/{1}{2}.mrc'.format(
                    corrected_path, file_name, options.sum_suffix
                    )
            frames_name = '{0}/{1}{2}.mrc'.format(
                    corrected_path, file_name, options.frames_suffix
                    )
            micrograph_name_skip = '{0}/{1}{2}.mrc'.format(
                    uncorrected_path, file_name, options.sum_suffix
                    )
            frames_name_skip = '{0}/{1}{2}.mrc'.format(
                    uncorrected_path, file_name, options.frames_suffix
                    )
        shift_name = '{0}/{1}{2}.txt'.format(
                shift_path, file_name, options.shift_suffix
                )
        frc_name = '{0}/{1}{2}.txt'.format(
                frc_path, file_name, options.frc_suffix
                )
        if not options.summovie_ready:
            temp_name = '{0}/{1}{2}.mrc'.format(
                    temp_path, file_name, options.sum_suffix
                    )
        else:
            temp_name = inputfile
        log_name = '{0}/{1}.log'.format(
                log_path, file_name
                )
        error_name = '{0}/{1}.err'.format(
                log_path, file_name
                )
        # Append the names to the lists
        micrograph_list.append('{0}{1}.mrc'.format(file_name, options.sum_suffix))
        shift_list.append(shift_name)
        if options.save_frames:
            frames_list.append('{0}{1}.mrc'.format(file_name, options.frames_suffix))


        # First build the summovie command
        if not options.skip_dose_filter:
            summovie_command = create_summovie_command(
                temp_name,
                micrograph_name,
                shift_name,
                frames_name,
                options
                )
            summovie_command_skip = create_summovie_command(
                temp_name,
                micrograph_name_skip,
                shift_name,
                frames_name_skip,
                options,
                skip=True
                )
        else:
            summovie_command = create_summovie_command(
                temp_name,
                micrograph_name,
                shift_name,
                frames_name,
                options
                )


        # Export the number of threads
        export_threads_command = []

        # Export
        export_threads_command.append('export')
        # Nr of threads
        export_threads_command.append('OMP_NUM_THREADS={0}'.format(
            options.nr_threads
            ))

        if not options.summovie_ready:
            # Do a e2proc3d.py
            e2proc3d_command = []

            # e2proc3d
            e2proc3d_command.append('e2proc3d.py')
            # inputfile
            e2proc3d_command.append('{0}'.format(inputfile))
            # outputfile
            e2proc3d_command.append('{0}'.format(temp_name))


        # Translate the command to single strings
        if not options.summovie_ready:
            e2proc3d_command = r' '.join(e2proc3d_command)
        export_threads_command = r' '.join(export_threads_command)
        summovie_command = '\n'.join(summovie_command)
        if not options.skip_dose_filter:
            summovie_command_skip = '\n'.join(summovie_command_skip)

        # Build full command
        if not options.summovie_ready:
            if not options.skip_dose_filter:
                full_command = r'{0}; {1}; echo "{2}" | {3}'.format(
                        export_threads_command,
                        e2proc3d_command,
                        summovie_command,
                        summovie_path
                        )
                full_command_skip = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command_skip,
                        summovie_path
                        )
            else:
                full_command = r'{0}; {1}; echo "{2}" | {3}'.format(
                        export_threads_command,
                        e2proc3d_command,
                        summovie_command,
                        summovie_path
                        )
        else:
            if not options.skip_dose_filter:
                full_command = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command,
                        summovie_path
                        )
                full_command_skip = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command_skip,
                        summovie_path
                        )
            else:
                full_command = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command,
                        summovie_path
                        )

        # Remove temp summovie files
        temp_summovie_files = glob('.UnBlur*')
        for entry in temp_summovie_files:
            remove(entry)

        with open(log_name, 'w') as f:
            with open(error_name, 'w') as e:
                # Execute Command
                if not options.skip_dose_filter:
                    subprocess.Popen(
                        [full_command], shell=True,
                        stdout=f,
                        stderr=e
                        ).wait()

                    # Remove temp summovie files
                    temp_summovie_files = glob('.UnBlur*')
                    for entry in temp_summovie_files:
                        remove(entry)

                    subprocess.Popen(
                        [full_command_skip], shell=True,
                        stdout=f,
                        stderr=e
                        ).wait()
                else:
                    subprocess.Popen(
                        [full_command], shell=True,
                        stdout=f,
                        stderr=e
                        ).wait()

        # Remove temp summovie files
        temp_summovie_files = glob('.UnBlur*')
        for entry in temp_summovie_files:
            remove(entry)
        if not options.summovie_ready:
            if path.exists(temp_name):
                # Remove temp file
                remove(temp_name)
            else:
                print('Error with file:\n{0}'.format(inputfile))

        time_list.append(time.time() - t1)

        # Do progress output
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
        print(
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
            )


    # Write micrograph and shift list
    with open('{0}/summovie_micrographs.txt'.format(output_dir), 'w') as f:
        for entry in sorted(micrograph_list):
            f.write('{0}\n'.format(entry))

    with open('{0}/summovie_shiftfiles.txt'.format(output_dir), 'w') as f:
        for entry in sorted(shift_list):
            f.write('{0}\n'.format(entry))

    if options.save_frames:
        with open('{0}/summovie_frames.txt'.format(output_dir), 'w') as f:
            for entry in sorted(frames_list):
                f.write('{0}\n'.format(entry))


def create_summovie_command(
        temp_name,
        micrograph_name,
        shift_name,
        frames_name,
        options,
        skip=False
        ):

    # Command list
    summovie_command = []

    # Input file
    summovie_command.append('{0}'.format(temp_name))
    # Number of frames
    summovie_command.append('{0}'.format(options.nr_frames))
    # Sum file
    summovie_command.append(micrograph_name)
    # Shift file
    summovie_command.append(shift_name)
    # Pixel size
    summovie_command.append('{0}'.format(options.pixel_size))
    # Dose correction
    if options.skip_dose_filter or skip:
        summovie_command.append('NO')
    else:
        summovie_command.append('YES')
        # Exposure per frame
        summovie_command.append('{0}'.format(options.exposure_per_frame))
        # Acceleration voltage
        summovie_command.append('{0}'.format(options.voltage))
        # Pre exposure
        summovie_command.append('{0}'.format(options.pre_exposure))
    # Save frames
    if not options.save_frames:
        summovie_command.append('NO')
    else:
        summovie_command.append('YES')
        # Frames output
        summovie_command.append('{0}'.format(frames_name))
    # Expert mode
    if not options.expert_mode:
        summovie_command.append('NO')
    else:
        summovie_command.append('YES')
        # FRC file
        summovie_command.append('{0}'.format(frc_name))
        # Minimum shift for initial search
        summovie_command.append('{0}'.format(options.shift_initial))
        # Outer radius shift limit
        summovie_command.append('{0}'.format(options.shift_radius))
        # B-Factor to Apply
        summovie_command.append('{0}'.format(options.b_factor))
        # Half-width vertical
        summovie_command.append('{0}'.format(options.fourier_vertical))
        # Half-width horizontal
        summovie_command.append('{0}'.format(options.fourier_horizontal))
        # Termination shift threshold
        summovie_command.append('{0}'.format(options.shift_threshold))
        # Maximum iterations
        summovie_command.append('{0}'.format(options.iterations))
        # Restore noise power
        if options.dont_restore_noise:
            summovie_command.append('NO')
        else:
            summovie_command.append('YES')
        # Verbose output
        if options.verbose:
            summovie_command.append('YES')
        else:
            summovie_command.append('NO')

    return summovie_command


if __name__ == '__main__':
    main()
