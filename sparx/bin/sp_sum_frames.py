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

from sparx import EMData, EMUtil, Region, fshift, get_im, read_text_row, fft, Util
from sys import argv
from os import path
from numpy import genfromtxt
import numpy
from global_def import SPARXVERSION, ERROR
from optparse import OptionParser
import global_def
import time

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

    unblur_image = EMData('unblur_small.mrc')
    unblur2_image = EMData('unblur2.mrc')
    t1 = time.time()
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

#    ########## Translation in fourier space using fshift
#    if not options.skip_alignment:
#        # Import shift file in angstrom by row 
#        sx, sy = genfromtxt(options.shift_file)
#
#        # Translate angstrom to pixel
#        sx /= options.pixel_size
#        sy /= options.pixel_size
#
#        # Translate the frames and sum them
#        for i, x, y in zip(
#                range(first, last), 
#                sx[first:last], 
#                sy[first:last]
#                ):
#            # Get the single frame
#            frame = input_stack.get_clip(Region(0, 0, i, nx, ny, 1))
#            frame = frame.do_fft()
#            # Transform the frame into fourier space
#            frame = fshift(frame, x, y)
#            # Add to a sum image
#            if i == first:
#                output_stack = frame
#            else:
#                output_stack.add(frame)
#        output_stack = output_stack.do_ift()
#    else:
#        # Translate the frames and sum them
#        for i in range(options.first, last): 
#            # Get the single frame
#            frame = input_stack.get_clip(Region(0, 0, i, round(nx), round(ny), 1))
#            # Add to a sum image
#            output_stack = output_stack + frame
#
#    # Write output
#    output_stack.write_image(output_name)
#    print('fshift:', time.time()-t1)
#
#
#    diff = output_stack - unblur_image
#    diff.write_image('diff.mrc')
#
#    #### Using pawels method
#    t1 = time.time()
#
#    movie = get_im(input_name)
#    nx = movie.get_xsize()
#    ny = movie.get_ysize()
#    nz = movie.get_zsize()
#
#    shifts = read_text_row(options.shift_file)
#    if len(shifts[-2]) != nz:
#        print('nope')
#
#    o = EMData(nx, ny, 1, False)
#    for i in xrange(nz):
#        t = fft(movie.get_clip( Region( 0, 0, i, nx, ny, 1) ))
#        Util.add_img(o, fshift(t, shifts[-2][i]/options.pixel_size, shifts[-1][i]/options.pixel_size))
#
#    fft(o).write_image('pawels_suggestion.mrc')
#    # Write output
#    print('fshift pawel:', time.time()-t1)
#    diff = fft(o) - unblur_image
#    diff.write_image('diff_pawel.mrc')


    ### DO IT ****** MANUAL!!!!
    t1 = time.time()

    movie = get_im(input_name)
    nx = movie.get_xsize()
    ny = movie.get_ysize()
    nz = movie.get_zsize()
    print(nx, ny, nz)
    
    shifts = read_text_row(options.shift_file)
    sx = shifts[-2]
    sy = shifts[-1]
   
    o = EMData(nx, ny, 1, False)
    for i in xrange(nz):
        t = fft(movie.get_clip( Region(0, 0, i, nx, ny, 1) ))
        f = fourier_shifterino(t, sx[i], sy[i], nx, ny)
        Util.add_img(o, f)

    fft(o).write_image('own_stuff.mrc')
    # Write output
    print('own: ', time.time() - t1)
    diff = fft(o) - unblur_image
    diff.write_image('diff_own.mrc')




    global_def.BATCH = False


def fourier_shifterino(frame, xshift, yshift, nx, ny):
    zshift=0
    nz = frame.get_zsize()
    ny = frame.get_ysize()
    nx = frame.get_xsize()
    if nz > 1:
        nzp = nz
    else:
        nzp = 1
    if ny > 1:
        nyp = ny
    else:
        nyp = 1
    nxp = nx
    lsd2 = (nxp + 2 - nxp%2) / 2
    nzp2 = nzp/2
    nyp2 = nyp/2
    print(frame.get_complex_at(0, 0, 0))
    print(frame.get_data())
#    for iz in xrange(1, nzp+1):
#        jz=iz-1 
#        if jz>nzp2:
#            jz=jz-nzp
#        for iy in xrange(1, nyp+1):
#            jy=iy-1
#            if jy>nyp2:
#                jy=jy-nyp
#            for ix in xrange(1, lsd2+1):
#                jx=ix-1
#                b = frame.get_complex_at(ix, iy, iz)
#                a = b * numpy.exp(
#                    -2*float(numpy.pi)*1j*(
#                        xshift*jx/nx + yshift*jy/ny+ zshift*jz/nz
#                        )
#                    )
#                frame.set_complex_at(ix, iy, iz, a)

    return frame

if __name__ == '__main__':
    main()
