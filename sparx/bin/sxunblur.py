#! /usr/bin/env python
#
# Copyright (C) 2016  Markus Stabrin (markus.stabrin@mpi-dortmund.mpg.de)
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

from sys import argv
from os import mkdir, path, system
from glob import glob
import global_def
from global_def import SPARXVERSION, ERROR
from optparse import OptionParser, SUPPRESS_HELP


def main():

    # Parse the Options
    progname = path.basename(argv[0])
    usage = progname + """ unblur input_image output
    --nr_frames=nr_frames
    --pixel_size=pixel_size
    --dose_filter
    --exposure_per_frame=exposure_per_frame
    --voltage=voltage
    --pre_exposure=pre_exposure
    --save_frames --expert_mode
    --shift_initial=shift_initial
    --shift_radius=shift_radius
    --b_factor=b_factor
    --fourier_vertical=fourier_vertical
    --fourier_horizontal=fourier_horizontal
    --shift_threshold=shift_threshold
    --iterations=iterations
    --restore_noise
    --verbose
    --filter_sum
    --lowpass=lowpass
    --highpass=highpass
    --remove_sum'

    sxunblur exists in non-MPI version.

    Just shift data.

    sxunblur.py directory_to_unblur directory/prefix*suffix.mrc output_directory
    --nr_frames=25 --pixel_size=1.19 --remove_sum

    Shift data with aligned sum files, filtered sum files and aligned frames.

    sxunblur.py directory_to_unblur directory/prefix*suffix.mrc output_directory
    --nr_frames=25 --pixel_size=1.19 --save_frames --filter_sum
    --lowpass=0.033 --highpass=0.00033

    Dose filter and Expert Options

    sxunblur.py directory_to_unblur directory/prefix*suffix.mrc output_directory
    --nr_frames=25 --pixel_size=1.19 --dose_filter --exposure_per_frame=1.0
    --voltage=300.0 --pre_exposure=0.0 --save_frames --expert_mode
    --shift_initial=2.0 --shift_radius=200.0 --b_factor=1500.0
    --fourier_vertical=1 --fourier_horizontal=1 --shift_threshold=0.1
    --iterations=10 --restore_noise --verbose --filter_sum --lowpass=0.033
    --highpass=0.00033
    """

    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option('--nr_frames',                   type='int',                  default=3,                                  help='number of frames in the set of micrographs')
    parser.add_option('--sum_suffix',                  type='str',                  default='_sum',                           help=SUPPRESS_HELP)
    parser.add_option('--shift_suffix',                  type='str',                    default='_shift',                        help=SUPPRESS_HELP)
    parser.add_option('--pixel_size',                    type='float',                  default=-1.0,                          help='pixel size [A]')
    parser.add_option('--dose_filter',                   action='store_true',       default=False,                        help='apply dose filter options')
    parser.add_option('--exposure_per_frame',    type='float',                  default=2.0,                           help='exposure per frame [e/A^2]')
    parser.add_option('--voltage',                         type='float',                    default=300.0,                    help='accelerate voltage [kV]')
    parser.add_option('--pre_exposure',                type='float',                    default=0.0,                       help='pre exposure amount [e/A^2]')
    parser.add_option('--save_frames',                  action='store_true',        default=False,                     help='save aligned frames')
    parser.add_option('--frames_suffix',                type='string',                  default='_sum_frames',       help=SUPPRESS_HELP)
    parser.add_option('--expert_mode',                 action='store_true',        default=False,                    help='set expert mode settings')
    parser.add_option('--frc_suffix',                       type='string',                  default='_frc',                     help=SUPPRESS_HELP)
    parser.add_option('--shift_initial',                    type='float',                   default=2.0,                        help='minimum shift for inital search [A]')
    parser.add_option('--shift_radius',                    type='float',                    default=200.0,                   help='outer radius shift limit [A]')
    parser.add_option('--b_factor',                         type='float',                   default=1500.0,                  help='b-factor to appy to image [A^2]')
    parser.add_option('--fourier_vertical',              type='int',                    default=1,                             help='half-width of central vertical line of fourier mask')
    parser.add_option('--fourier_horizontal',         type='int',                   default=1,                               help='half-width of central horizontal line of fourier mask')
    parser.add_option('--shift_threshold',               type='float',                  default=0.1,                         help='termination shift threshold')
    parser.add_option('--iterations',                       type='int',                     default=10,                           help='maximum number of iterations')
    parser.add_option('--restore_noise',                action='store_true',         default=False,                     help='restore noise power')
    parser.add_option('--verbose',                          action='store_true',         default=False,                    help='verbose output')
    parser.add_option('--filter_sum',                       action='store_true',        default=False,                     help='filter the output images')
    parser.add_option('--lowpass',                          type='float',                   default=0.033,                     help='apply a lowpass filter: abolute frequency')
    parser.add_option('--highpass',                         type='float',                   default=0.00033,                help='apply a highpass filter: abolute frequency')
    parser.add_option('--remove_sum',                   action='store_true',        default=False,                     help='remove the calculated sum files')

    # list of the options and the arguments
    (options, args) = parser.parse_args(argv[1:])

    global_def.BATCH = True

    # If there arent enough arguments, stop the script
    if len(args) != 3:
        ERROR("see usage " + usage, 1)

    # Convert the realtive parts to absolute ones
    unblur_path = path.realpath(args[0])
    input_image = path.realpath(args[1])
    output_dir = path.realpath(args[2])

    # If the unblur executable file does not exists, stop the script
    if not path.exists(unblur_path):
        ERROR(
            'Unblur directory does not exist, please change' +
            ' the name and restart the program.', 1
            )

    # If the output directory exists, stop the script
    if path.exists(output_dir):
        ERROR(
            'Output directory exists, please change' +
            ' the name and restart the program.', 1
            )

    # If the input file does not exists, stop the script
    fileList = glob(input_image)

    if not fileList:
        ERROR(
            'Input file does not exist, please change' +
            ' the name and restart the program.', 1
            )

    # Split the path of the image name at the "/" Characters.
    # The last entry contains the micrograph name.
    # Split the micrograph name at the wildcard character for the
    # prefix and suffix.
    input_split = input_image.split('/')
    input_name = input_split[-1].split('*')

    # Check, if there is an prefix and suffix.
    if len(input_name) == 2:
        input_suffix = input_name[1]
    else:
        if '*' in input_split[-1]:
            input_suffix = input_name[0]
        else:
            input_suffix = '.mrc'

    if len(input_split) != 1:
        input_dir = input_image[:-len(input_split[-1])]
    else:
        input_dir = ''

    # Create output directorys
    mkdir('{:s}'.format(output_dir))
    if not path.exists('{:s}/Micrographs'.format(output_dir)):
        mkdir('{:s}/Micrographs'.format(output_dir))
    if not path.exists('{:s}/Shift'.format(output_dir)):
        mkdir('{:s}/Shift'.format(output_dir))
    if not path.exists('{:s}/Filtered'.format(output_dir)) \
            and options.filter_sum:
        mkdir('{:s}/Filtered'.format(output_dir))
    if not path.exists('{:s}/FRC'.format(output_dir)) \
            and options.expert_mode:
        mkdir('{:s}/FRC'.format(output_dir))

    # Create sh script
    create_sh_script(
        unblur_path=unblur_path,
        input_image=input_image,
        input_dir=input_dir,
        output_dir=output_dir,
        input_suffix=input_suffix,
        options=options
        )

    # Start sh script
    system('sh {:s}/scriptUnblur.sh'.format(output_dir))

    global_def.BATCH = False


def create_sh_script(
            unblur_path, input_image, output_dir,
            input_dir, input_suffix, options
            ):
        """Create the sh script for starting unblur"""
        strSh = ''

        # To make sure it is a bash script
        strSh += '#!/bin/bash\n\n'

        # Create a file list of all files
        strSh += 'fileList=$(ls {:s})\n'.format(
            input_image
            )

        # Create folders
        strSh += 'mkdir {:s}/Micrographs\n'.format(output_dir)

        strSh += 'mkdir {:s}/Shift\n'.format(output_dir)

        if options.filter_sum:
            strSh += 'mkdir {:s}/Filtered\n'.format(output_dir)

        if options.expert_mode:
            strSh += 'mkdir {:s}/FRC\n'.format(output_dir)

        # Loop over all files
        strSh += '\nfor file in $fileList\ndo\n\n'

        strSh += 'baseName=${{file%{:s}}}\n'.format(input_suffix)
        strSh += 'baseName=${{baseName#{:s}}}\n'.format(input_dir)

        # Create a temporary file to work with to prevent format issues
        strSh += 'e2proc3d.py $file tempUnblur.mrc\n\n'

        # Remove some temporary files that unblur makes
        strSh += 'rm .UnBlur*\n'

        # Start Unblur
        strSh += '{:s} << eof\n'.format(unblur_path)

        # Input File
        strSh += 'tempUnblur.mrc\n'
        # Number of Frames
        strSh += '{:d}\n'.format(options.nr_frames)
        # Sum File
        strSh += '{:s}/Micrographs/${{baseName}}{:s}.mrc\n'.format(
            output_dir,
            options.sum_suffix
            )
        # Shift File
        strSh += '{:s}/Shift/${{baseName}}{:s}.txt\n'.format(
            output_dir,
            options.shift_suffix
            )
        # Pixel Size
        strSh += '{:f}\n'.format(options.pixel_size)

        if options.dose_filter:
            # Say yes to Dose Filtering
            strSh += 'YES\n'
            # Exposure per Frame
            strSh += '{:f}\n'.format(options.exposure_per_frame)
            # Acceleration Voltage
            strSh += '{:f}\n'.format(options.voltage)
            # Pre Exposure
            strSh += '{:f}\n'.format(options.pre_exposure)
        else:
            # Say no to Dose Filtering
            strSh += 'NO\n'

        if options.save_frames:
            # Say yes to Save Frames
            strSh += 'YES\n'
            # Frames file
            strSh += '{:s}/Micrographs/${{baseName}}{:s}.mrc\n'.format(
                output_dir,
                options.frames_suffix
                )
        else:
            # Say no to Save Frames
            strSh += 'NO\n'

        if options.expert_mode:
            # Say yes to Expert Mode
            strSh += 'YES\n'
            # FRC File
            strSh += '{:s}/FRC/${{baseName}}{:s}.txt\n'.format(
                output_dir,
                options.frc_suffix
                )
            # Minimum Shift for initial search
            strSh += '{:f}\n'.format(options.shift_initial)
            # Outer Radius Shift Limit
            strSh += '{:f}\n'.format(options.shift_radius)
            # B-Factor to Apply
            strSh += '{:f}\n'.format(options.b_factor)
            # Half-Width Vertical
            strSh += '{:d}\n'.format(options.fourier_vertical)
            # Hald-Width Horizontal
            strSh += '{:d}\n'.format(options.fourier_horizontal)
            # Termination Shift Threshold
            strSh += '{:f}\n'.format(options.shift_threshold)
            # Maximum Iterations
            strSh += '{:d}\n'.format(options.iterations)
            # Restore Noise Power
            if options.restore_noise:
                # Say yes to Restore Noise Power
                strSh += 'YES\n'
            else:
                # Say no to Restore Noise Power
                strSh += 'NO\n'
            # Verbose Output
            if options.verbose:
                # Say yes to Verbose Output
                strSh += 'YES\n'
            else:
                # Say no to Verbose Output
                strSh += 'NO\n'
        else:
            # Say no to Expert Mode
            strSh += 'NO\n'

        # Enf of file reached
        strSh += 'eof\n\n'

        # Remove temporary file
        strSh += 'rm tempUnblur.mrc\n'

        # Remove some temporary files that unblur makes
        strSh += 'rm .UnBlur*\n\n'

        if options.filter_sum:
            # Filter Images
            lowpass_angstrom = options.pixel_size / options.lowpass
            highpass_angstrom = options.pixel_size / options.highpass
            strSh += \
                'e2proc3d.py {:s}/Micrographs/${{baseName}}{:s}.mrc '.format(
                    output_dir,
                    options.sum_suffix
                    )
            strSh += '{:s}/Filtered/${{baseName}}_lp{:d}_hp{:d}.mrc ' \
                .format(
                    output_dir,
                    int(lowpass_angstrom),
                    int(highpass_angstrom)
                    )
            strSh += '--process=filter.lowpass.gauss:cutoff_freq={:f} '.format(
                options.lowpass
                )
            strSh += '--process=filter.highpass.gauss:cutoff_freq={:f}\n\n' \
                .format(
                    options.highpass
                    )

        if options.remove_sum:
            # Remove sum files
            strSh += 'rm {:s}/Micrographs/${{baseName}}{:s}.mrc\n'.format(
                output_dir,
                options.sum_suffix
                )

        # Done
        strSh += 'done\n'

        # Write Output
        with open('{:s}/scriptUnblur.sh'.format(output_dir), 'w') as f:
            f.write(strSh)

if __name__ == '__main__':
    main()
