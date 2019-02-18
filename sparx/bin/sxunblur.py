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
    usage = progname + """ unblur_path input_micrograph_pattern output_directory
    --summovie_path
    --selection_list
    --nr_frames=nr_frames
    --pixel_size=pixel_size
    --voltage=voltage
    --exposure_per_frame=exposure_per_frame
    --pre_exposure=pre_exposure
    --nr_threads
    --save_frames
    --skip_dose_filter
    --expert_mode
    --shift_initial=shift_initial
    --shift_radius=shift_radius
    --b_factor=b_factor
    --fourier_vertical=fourier_vertical
    --fourier_horizontal=fourier_horizontal
    --shift_threshold=shift_threshold
    --iterations=iterations
    --dont_restore_noise
    --verbose

    sxunblur exists only in non-MPI version.

    Perform unblur and with dose filtering and summovie without dose filtering.

    sxunblur.py ~/my_app/unblur 'movies/micrograph_*_frames.mrc' outdir_unblur
    --summovie_path=~/my_app/summovie
    --nr_frames=25 --pixel_size=1.19 --exposure_per_frame=1.0
    --voltage=300.0 --pre_exposure=0.0 --nr_threads=1

    Perform unblur with dose filtering and summovie without dose filtering with selection list.

    sxunblur.py ~/my_app/unblur 'movies/micrograph_*_frames.mrc' outdir_unblur
    --summovie_path=~/my_app/summovie
    --selection_list=selected_micrograph_file
    --nr_frames=25 --pixel_size=1.19 --exposure_per_frame=1.0
    --voltage=300.0 --pre_exposure=0.0 --nr_threads=1

    Perform unblur without dose filtering.

    sxunblur.py ~/my_app/unblur 'movies/micrograph_*_frames.mrc' outdir_unblur
    --nr_frames=25 --pixel_size=1.19 --skip_dose_filter --nr_threads=1

    Perform unblur without dose filtering and save the frames.

    sxunblur.py ~/my_app/unblur 'movies/micrograph_*_frames.mrc' outdir_unblur
    --nr_frames=25 --pixel_size=1.19 --skip_dose_filter --save_frames --nr_threads=1

    Perform unblur with dose filtering and summovie without dose filtering with all options.

    sxunblur.py ~/my_app/unblur 'movies/micrograph_*_frames.mrc' outdir_unblur
    --summovie_path=~/my_app/summovie
    --nr_frames=25 --pixel_size=1.19 --exposure_per_frame=1.0
    --voltage=300.0 --pre_exposure=0.0 --save_frames --expert_mode
    --shift_initial=2.0 --shift_radius=200.0 --b_factor=1500.0
    --fourier_vertical=1 --fourier_horizontal=1 --shift_threshold=0.1
    --iterations=10 --verbose --nr_threads=1
    """

    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option('--summovie_path',      type='str',          default='',        help='summovie executable path (SPHIRE specific): Specify the file path of summovie executable. (default none)')
    parser.add_option('--selection_list',     type='str',          default='',        help='Micrograph selecting list (SPHIRE specific): Specify a name of micrograph selection list text file. The file extension must be \'.txt\'. If this is not provided, all files matched with the micrograph name pattern will be processed. (default none)')
    parser.add_option('--nr_frames',          type='int',          default=3,         help='Number of movie frames: The number of movie frames in each input micrograph. (default 3)')
    parser.add_option('--sum_suffix',         type='str',          default='_sum',    help=SUPPRESS_HELP)
    parser.add_option('--shift_suffix',       type='str',          default='_shift',  help=SUPPRESS_HELP)
    parser.add_option('--pixel_size',         type='float',        default=-1.0,      help='Pixel size [A]: The pixel size of input micrographs. (default required float)')
    parser.add_option('--voltage',            type='float',        default=300.0,     help='Microscope voltage [kV]: The acceleration voltage of microscope used for imaging. (default 300.0)')
    parser.add_option('--exposure_per_frame', type='float',        default=2.0,       help='Per frame exposure [e/A^2]: The electron dose per frame in e/A^2. (default 2.0)')
    parser.add_option('--pre_exposure',       type='float',        default=0.0,       help='Pre-exposure [e/A^2]: The electron does in e/A^2 used for exposure prior to imaging .(default 0.0)')
    parser.add_option('--nr_threads',         type='int',          default=1,         help='Number of threads: The number of threads unblur can use. The higher the faster, but it requires larger memory. (default 1)')
    parser.add_option('--save_frames',        action='store_true', default=False,     help='Save aligned movie frames: Save aligned movie frames. This option slows down the process. (default False)')
    parser.add_option('--frames_suffix',      type='string',       default='_frames', help=SUPPRESS_HELP)
    parser.add_option('--skip_dose_filter',   action='store_true', default=False,     help='Skip dose filter step: With this option, voltage, exposure per frame, and pre exposure will be ignored. (default False)')
    parser.add_option('--expert_mode',        action='store_true', default=False,     help='Use expert mode: Requires initial shift, shift radius, b-factor, fourier_vertical, fourier_horizontal, shift threshold, iterations, restore noise, and verbosity options. (default False)')
    parser.add_option('--frc_suffix',         type='string',       default='_frc',    help=SUPPRESS_HELP)
    parser.add_option('--shift_initial',      type='float',        default=2.0,       help='Minimum shift for initial search [A] (expert mode): Effective with unblur expert mode. (default 2.0)')
    parser.add_option('--shift_radius',       type='float',        default=200.0,     help='Outer radius shift limit [A] (expert mode): Effective with unblur expert mode. (default 200.0)')
    parser.add_option('--b_factor',           type='float',        default=1500.0,    help='Apply B-factor to images [A^2] (expert mode): Effective with unblur expert mode. (default 1500.0)')
    parser.add_option('--fourier_vertical',   type='int',          default=1,         help='Vertical Fourier central mask size (expert mode): The half-width of central vertical line of Fourier mask. Effective with unblur expert mode. (default 1)')
    parser.add_option('--fourier_horizontal', type='int',          default=1,         help='Horizontal Fourier central mask size (expert mode): The half-width of central horizontal line of Fourier mask. Effective with unblur expert mode. (default 1)')
    parser.add_option('--shift_threshold',    type='float',        default=0.1,       help='Termination shift threshold (expert mode): Effective with unblur expert mode. (default 0.1)')
    parser.add_option('--iterations',         type='int',          default=10,        help='Maximum iterations (expert mode): Effective with unblur expert mode. (default 10)')
    parser.add_option('--dont_restore_noise', action='store_true', default=False,     help='Do not restore noise power (expert mode): Effective with unblur expert mode. (default False)')
    parser.add_option('--verbose',            action='store_true', default=False,     help='Verbose (expert mode): Effective with unblur expert mode. (default False)')
    parser.add_option('--unblur_ready',       action='store_true', default=False,     help=SUPPRESS_HELP)

    # list of the options and the arguments
    (options, args) = parser.parse_args(argv[1:])

    global_def.BATCH = True

    # If there arent enough arguments, stop the script
    if len(args) != 3:
        ERROR("see usage " + usage, 1)

    # Convert the realtive parts to absolute ones
    unblur_path = path.realpath(args[0]) # unblur_path
    input_image = path.realpath(args[1]) # input_micrograph_pattern
    output_dir = path.realpath(args[2])  # output_directory

    # If the unblur executable file does not exists, stop the script
    if not path.exists(unblur_path):
        ERROR(
            'Unblur directory does not exist, please change' +
            ' the name and restart the program.', 'sxunblur.py', 1
            )

    # If the output directory exists, stop the script
    if path.exists(output_dir):
        ERROR(
            'Output directory exists, please change' +
            ' the name and restart the program.', 'sxunblur.py', 1
            )

    # If the input file does not exists, stop the script
    file_list = glob(input_image)

    if not file_list:
        ERROR(
            'Input file does not exist, please change' +
            ' the name and restart the program.', 'sxunblur.py', 1
            )

    # If the skip_dose_filter option is false, the summovie path is necessary
    if not options.skip_dose_filter and not path.exists(options.summovie_path):
        ERROR(
            'Path to the SumMovie executable is necessary when dose weighting is performed.',
            'sxunblur.py', 1
            )


    # Output paths
    corrected_path = '{:s}/corrsum_dose_filtered'.format(output_dir)
    uncorrected_path = '{:s}/corrsum'.format(output_dir)
    shift_path = '{:s}/shift'.format(output_dir)
    frc_path = '{:s}/frc'.format(output_dir)
    log_path = '{:s}/logfiles'.format(output_dir)
    temp_path = '{0}/temp'.format(output_dir)

    # Split the path of the image name at the "/" Characters.
    # The last entry contains the micrograph name.
    # Split the micrograph name at the wildcard character for the
    # prefix and suffix.
    input_split = input_image.split('/')
    input_name = input_split[-1].split('*')

    # Get the input directory
    if len(input_split) != 1:
        input_dir = input_image[:-len(input_split[-1])]
    else:
        input_dir = ''

    # Create output directorys
    if not path.exists(output_dir):
        mkdir(output_dir)
    if not path.exists(uncorrected_path):
        mkdir(uncorrected_path)
    if not path.exists(shift_path):
        mkdir(shift_path)
    if not path.exists(corrected_path):
        if not options.skip_dose_filter:
            mkdir(corrected_path)
    if not path.exists(frc_path):
        if options.expert_mode or not options.skip_dose_filter:
            mkdir(frc_path)
    if not path.exists(temp_path) and not options.unblur_ready:
        mkdir(temp_path)
    if not path.exists(log_path):
        mkdir(log_path)

    # Run unblur
    run_unblur(
        unblur_path=unblur_path,
        input_image=input_image,
        input_dir=input_dir,
        output_dir=output_dir,
        corrected_path=corrected_path,
        uncorrected_path=uncorrected_path,
        shift_path=shift_path,
        frc_path=frc_path,
        temp_path=temp_path,
        log_path=log_path,
        file_list=file_list,
        options=options
        )

    if not options.unblur_ready:
        # Remove temp folder
        for entry in glob('{0}/*'.format(temp_path)):
            remove(entry)
        rmdir(temp_path)

    print('All Done!')

    global_def.BATCH = False


def run_unblur(
        unblur_path,
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
    mic_list = options.selection_list
    if mic_list:
        # Import list file
        try:
            set_selection = genfromtxt(mic_list, dtype=None)
        except TypeError:
            ERROR('no entrys in list file {0}'.format(mic_list), 'sxunblur.py', 1)
        # List of files which are in pattern and list
        file_list = [
                entry for entry in file_list \
                if entry[len(input_dir):] in set_selection and \
                path.exists(entry)
                ]
        # If no match is there abort
        if len(file_list) == 0:
            ERROR(
                'no files in {0} matched the file pattern:\n'.format(mic_list), 'sxunblur.py',
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
            frc_name = '{0}/{1}{2}.txt'.format(
                frc_path, file_name, options.frc_suffix
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
            frc_name = '{0}/{1}{2}.txt'.format(
                frc_path, file_name, options.frc_suffix
                )
            frc_summovie_name = '{0}/{1}_summovie{2}.txt'.format(
                frc_path, file_name, options.frc_suffix
                )
        shift_name = '{0}/{1}{2}.txt'.format(
                shift_path, file_name, options.shift_suffix
                )
        if not options.unblur_ready:
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


        # First build the unblur/summovie command
        if not options.skip_dose_filter:
            unblur_command = create_unblur_command(
                temp_name,
                micrograph_name,
                shift_name,
                frames_name,
                frc_name,
                options
                )
            # Options for the summovie command
            options_summovie = {
                'first': 1,
                'last': -1,
                'nr_frames': options.nr_frames,
                'pixel_size': options.pixel_size,
                'exposure_per_frame': options.exposure_per_frame,
                'voltage': options.voltage,
                'pre_exposure': options.pre_exposure,
                'dont_restore_noise': options.dont_restore_noise,
                'apply_dose_filter': False
            }
            summovie_command = create_summovie_command(
                temp_name,
                micrograph_name_skip,
                shift_name,
                frc_summovie_name,
                options_summovie
                )
        else:
            unblur_command = create_unblur_command(
                temp_name,
                micrograph_name,
                shift_name,
                frc_name,
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

        if not options.unblur_ready:
            # Do a e2proc3d.py
            e2proc3d_command = []

            # e2proc3d
            e2proc3d_command.append('e2proc3d.py')
            # inputfile
            e2proc3d_command.append('{0}'.format(inputfile))
            # outputfile
            e2proc3d_command.append('{0}'.format(temp_name))


        # Translate the command to single strings
        if not options.unblur_ready:
            e2proc3d_command = r' '.join(e2proc3d_command)
        export_threads_command = r' '.join(export_threads_command)
        unblur_command = '\n'.join(unblur_command)
        if not options.skip_dose_filter:
            summovie_command = '\n'.join(summovie_command)

        # Build full command
        if not options.unblur_ready:
            if not options.skip_dose_filter:
                full_command = r'{0}; {1}; echo "{2}" | {3}'.format(
                        export_threads_command,
                        e2proc3d_command,
                        unblur_command,
                        unblur_path
                        )
                full_command_summovie = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command,
                        options.summovie_path
                        )
            else:
                full_command = r'{0}; {1}; echo "{2}" | {3}'.format(
                        export_threads_command,
                        e2proc3d_command,
                        unblur_command,
                        unblur_path
                        )
        else:
            if not options.skip_dose_filter:
                full_command = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        unblur_command,
                        unblur_path
                        )
                full_command_summovie = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        summovie_command,
                        options.summovie_path
                        )
            else:
                full_command = r'{0}; echo "{1}" | {2}'.format(
                        export_threads_command,
                        unblur_command,
                        unblur_path
                        )

        # Remove temp unblur files
        temp_unblur_files = glob('.UnBlur*')
        for entry in temp_unblur_files:
            remove(entry)

        # Remove temp summovie files
        temp_summovie_files = glob('.SumMovie*')
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

                    # Remove temp unblur files
                    temp_unblur_files = glob('.UnBlur*')
                    for entry in temp_unblur_files:
                        remove(entry)

                    # Remove temp summovie files
                    temp_summovie_files = glob('.SumMovie*')
                    for entry in temp_summovie_files:
                        remove(entry)

                    subprocess.Popen(
                        [full_command_summovie], shell=True,
                        stdout=f,
                        stderr=e
                        ).wait()
                else:
                    subprocess.Popen(
                        [full_command], shell=True,
                        stdout=f,
                        stderr=e
                        ).wait()

        # Remove temp unblur files
        temp_unblur_files = glob('.UnBlur*')
        for entry in temp_unblur_files:
            remove(entry)
        # Remove temp summovie files
        temp_summovie_files = glob('.SumMovie*')
        for entry in temp_summovie_files:
            remove(entry)

        if not options.unblur_ready:
            if path.exists(temp_name):
                # Remove temp file
                remove(temp_name)
            else:
                print(('Error with file:\n{0}'.format(inputfile)))

        # Check if SumMovie and UnBlur finished cleanly
        with open(log_name, 'r') as r:
            clean_summovie = False
            clean_unblur = False
            for line in r:
                if 'SumMovie finished cleanly.' in line:
                    clean_summovie = True
                if 'UnBlur finished cleanly.' in line:
                    clean_unblur = True

        if clean_unblur:
            print('UnBlur finished cleanly.')
        else:
            ERROR(
                'unblur error. check the logfile for more information: {0}'.format(
                    log_name
                    ), 'sxunblur.py', 0
                )

        if clean_summovie:
            print('SumMovie finished cleanly.')
        else:
            ERROR(
                'summovie error. check the logfile for more information: {0}'.format(
                    log_name
                    ), 'sxunblur.py', 0
                )

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


    # Write micrograph and shift list
    with open('{0}/unblur_micrographs.txt'.format(output_dir), 'w') as f:
        for entry in sorted(micrograph_list):
            f.write('{0}\n'.format(entry))

    with open('{0}/unblur_shiftfiles.txt'.format(output_dir), 'w') as f:
        for entry in sorted(shift_list):
            f.write('{0}\n'.format(entry))

    if options.save_frames:
        with open('{0}/unblur_frames.txt'.format(output_dir), 'w') as f:
            for entry in sorted(frames_list):
                f.write('{0}\n'.format(entry))


def create_unblur_command(
        temp_name,
        micrograph_name,
        shift_name,
        frames_name,
        frc_name,
        options
        ):

    # Command list
    unblur_command = []

    # Input file
    unblur_command.append('{0}'.format(temp_name))
    # Number of frames
    unblur_command.append('{0}'.format(options.nr_frames))
    # Sum file
    unblur_command.append(micrograph_name)
    # Shift file
    unblur_command.append(shift_name)
    # Pixel size
    unblur_command.append('{0}'.format(options.pixel_size))
    # Dose correction
    if options.skip_dose_filter:
        unblur_command.append('NO')
    else:
        unblur_command.append('YES')
        # Exposure per frame
        unblur_command.append('{0}'.format(options.exposure_per_frame))
        # Acceleration voltage
        unblur_command.append('{0}'.format(options.voltage))
        # Pre exposure
        unblur_command.append('{0}'.format(options.pre_exposure))
    # Save frames
    if not options.save_frames:
        unblur_command.append('NO')
    else:
        unblur_command.append('YES')
        # Frames output
        unblur_command.append('{0}'.format(frames_name))
    # Expert mode
    if not options.expert_mode:
        unblur_command.append('NO')
    else:
        unblur_command.append('YES')
        # FRC file
        unblur_command.append('{0}'.format(frc_name))
        # Minimum shift for initial search
        unblur_command.append('{0}'.format(options.shift_initial))
        # Outer radius shift limit
        unblur_command.append('{0}'.format(options.shift_radius))
        # B-Factor to Apply
        unblur_command.append('{0}'.format(options.b_factor))
        # Half-width vertical
        unblur_command.append('{0}'.format(options.fourier_vertical))
        # Half-width horizontal
        unblur_command.append('{0}'.format(options.fourier_horizontal))
        # Termination shift threshold
        unblur_command.append('{0}'.format(options.shift_threshold))
        # Maximum iterations
        unblur_command.append('{0}'.format(options.iterations))
        # Restore noise power
        if options.dont_restore_noise:
            unblur_command.append('NO')
        else:
            unblur_command.append('YES')
        # Verbose output
        if options.verbose:
            unblur_command.append('YES')
        else:
            unblur_command.append('NO')

    return unblur_command


if __name__ == '__main__':
    main()
