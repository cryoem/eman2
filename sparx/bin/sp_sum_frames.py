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

from sparx import EMData, EMUtil, Region
from sys import argv
from os import path
from numpy import genfromtxt
from global_def import SPARXVERSION, ERROR
from optparse import OptionParser
import global_def

def main():
    
    progname = path.basename(argv[0])
    usage = progname + """ input_stack output_stack
    --pixel_size
    --first
    --last
    --shift_file
    --skip_alignment

    sp_sum_frame exists in non-MPI version.
    """

    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option('--pixel_size', type=float, default=-1, help='Pixel size [A] (default -1)')
    parser.add_option('--first', type=int, default=0, help='First frame to use (default 0)')
    parser.add_option('--last', type=int, default=-1, help='Last frame to use (default -1)')
    parser.add_option('--shift_file', type=str, default='', help='Shift file for alignment (default "")')
    parser.add_option('--skip_alignment', action='store_true', default=False, help='Skip the alignment and just sum the frames (default False)')
    (options, args) = parser.parse_args(argv[1:])

    global_def.BATCH = True

    # If arguments are missing abort
    if len(args) != 2:
        ERROR("see usage:\n" + usage, 1)

    if options.pixel_size <= 0:
        ERROR("pixel size [A] needs to be greater then 0 (default -1)\n" + usage, 1)

    # Use realpath
    input_name = path.realpath(args[0])
    output_name = path.realpath(args[1])

    if not path.exists(input_name):
        ERROR("input image {0} does not exists".format(input_name), 1)

    if not path.exists(options.shift_file) and not \
            options.skip_alignment:
        ERROR(
            "given shift file {0} does not exists\n".format(
                options.shift_file
                ) + \
            usage, 
            1
            )

    sum_images(input_name, output_name, options)

    print('All Done!')

def sum_images(input_name, output_name, options):
    """Translate images and sum them"""
    
    #Import the input stack
    input_stack = EMData(input_name)

    # Get the dimensions of the input stack
    nx = input_stack['nx']
    ny = input_stack['ny']
    nz = input_stack['nz']

    # Check how many frames to use
    if nz < abs(options.last):
        ERROR(
                "Last option {0} is out of range: maximum |{1}|".format(
                options.last, nz-1
                ), 
            1)

    if nz < abs(options.first):
        ERROR(
                "First option {0} is out of range: maximum |{1}|".format(
                options.last, nz-1
                ), 
            1)

    # Get real indices
    if options.first < 0:
        first = nz + options.first
    else:
        first = options.first

    if options.last < 0:
        last = nz + options.last + 1
    else:
        last = options.last + 1

    if first >= last:
        ERROR("First option {0}={1} muss be smaller equals last option {2}={3}".format(
            options.first, first, options.last, last-1
            ),
            1)

    # Output volume 
    output_stack = EMData(nx, ny, 1)

    if not options.skip_alignment:
        # Import shift file in angstrom by row 
        sx, sy = genfromtxt(options.shift_file)

        # Translate angstrom to pixel
        sx /= options.pixel_size
        sy /= options.pixel_size

        # Translate the frames and sum them
        for i, x, y in zip(
                range(first, last), 
                sx[first:last], 
                sy[first:last]
                ):
            # Get the single frame
            frame = input_stack.get_clip(Region(0, 0, i, nx, ny, 1))
            # Transform the frame into fourier space
            frame = frame.do_fft()
            # If neccessary translate the images
            if x or y:
                frame.translate(x, y, 0.0)
            # Inverse fourier transform
            frame = frame.do_ift()
            # Add to a sum image
            output_stack = output_stack + frame

    else:
        # Translate the frames and sum them
        for i in range(options.first, last): 
            # Get the single frame
            frame = input_stack.get_clip(Region(0, 0, i, nx, ny, 1))
            # Add to a sum image
            output_stack = output_stack + frame

    # Write output
    output_stack.write_image(output_name)
    
    global_def.BATCH = False


if __name__ == '__main__':
    main()
