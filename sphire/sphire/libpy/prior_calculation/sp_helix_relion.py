"""
# Author: Markus Stabrin 2017-2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2017-2019 Max Planck Institute of Molecular Physiology

this software is issued under a joint bsd/gnu license. you may use the
source code in this file under either license. however, note that the
complete eman2 and sphire software packages have some gpl dependencies,
so you are responsible for compliance with the licenses of these packages
if you opt to use bsd licensing. the warranty disclaimer below holds
in either instance.

this complete copyright notice must be included in any revised version of the
source code. additional authorship citations may be added, but existing
author citations must be preserved.

this program is free software; you can redistribute it and/or modify
it under the terms of the gnu general public license as published by
the free software foundation; either version 2 of the license, or
(at your option) any later version.

this program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. see the
gnu general public license for more details.

you should have received a copy of the gnu general public license
along with this program; if not, write to the free software
foundation, inc., 59 temple place, suite 330, boston, ma  02111-1307 usa
"""

import numpy as np

def create_dtype_dict():
    """Create type dictionary"""
    # Open the dictionary file and convert it to a dictionary
    dict_as_string = \
        '_rlnCtfImage    |U200\n' + \
        '_rlnImageName   |U200\n' + \
        '_rlnMicrographName  |U200\n' + \
        '_rlnOriginalParticleName    |U200\n' + \
        '_rlnGroupName   |U200\n' + \
        '_rlnGroupNumber <i8\n' + \
        '_rlnClassNumber <i8\n' + \
        '_rlnNrOfSignificantSamples  <i8\n' + \
        '_rlnRandomSubset    <i8\n' + \
        '_rlnNrOfFrames  <i8\n' + \
        '_rlnHelicalTubeID  <i8\n' + \
        '_rlnVoltage <f8\n' + \
        '_rlnDefocusU    <f8\n' + \
        '_rlnDefocusV    <f8\n' + \
        '_rlnDefocusAngle    <f8\n' + \
        '_rlnSphericalAberration <f8\n' + \
        '_rlnDetectorPixelSize   <f8\n' + \
        '_rlnCtfFigureOfMerit    <f8\n' + \
        '_rlnMagnification   <f8\n' + \
        '_rlnAmplitudeContrast   <f8\n' + \
        '_rlnCoordinateX <f8\n' + \
        '_rlnCoordinateY <f8\n' + \
        '_rlnNormCorrection  <f8\n' + \
        '_rlnOriginX <f8\n' + \
        '_rlnOriginY <f8\n' + \
        '_rlnAngleRot    <f8\n' + \
        '_rlnAngleTilt   <f8\n' + \
        '_rlnAnglePsi    <f8\n' + \
        '_rlnAutopickFigureOfMerit   <f8\n' + \
        '_rlnLogLikeliContribution   <f8\n' + \
        '_rlnMaxValueProbDistribution    <f8\n' + \
        '_rlnParticleSelectZScore    <f8\n' + \
        '_rlnAngleRotPrior   <f8\n' + \
        '_rlnAngleTiltPrior  <f8\n' + \
        '_rlnAnglePsiPrior   <f8\n' + \
        '_rlnOriginXPrior    <f8\n' + \
        '_rlnOriginYPrior    <f8\n' + \
        '_rlnHelicalTrackLength  <f8\n'

    # Split the dictionary to a list
    string_list = dict_as_string.split()

    # Create a dictionary out of the list
    dtype_dict = {
        string_list[number]: string_list[number + 1]
        for number in range(0, len(string_list), 2)
        }

    return dtype_dict


def import_star_file(input_star_file):
    """Import the Data from the Star File as a structured Array"""
    # Load type dictionary
    dtype_dict = create_dtype_dict()

    # Create a list for the header information.  Use the built-in
    # function enumerate to go through the lines of the Star file.
    # Just save the names of the header information and stop after
    # the header is over.
    # If no header is found linenumber will stay False
    header_names = []
    linenumber = False
    for linenumber, line in enumerate(
            open(input_star_file, 'r')
            ):
        if line[0] == '_':
            header_names.append(line.split()[0])
        elif linenumber > 4:  # So data_ and loop_ won't abort the loop
            break
        elif linenumber > 50:
            break
        else:
            assert(linenumber <= 50)
            assert(linenumber <= 4)

    # Create a list for the dtype information.  Go through the
    # header_names list and append the dtype of the related column.
    # If there isn't an entry for the name yet it will be written to
    # the dictionary file and saved as float.
    # If no header_names are there set the column to ('column', '|U200')

    dtype_list = []
    if any(header_names):
        for names in header_names:
            try:
                dtype_list.append((names, dtype_dict[names]))
            except:
                dtype_list.append((names, '<f8'))
                dtype_dict.update({names: '<f8'})
    else:
        dtype_list.append(('column', '|U200'))
        linenumber = 0

    # Import the dataInput as a np structured Array.  Skip the lines
    # with header information and set the dtype_list.
    # If linenumber is False there is a input/output error.
    if linenumber:
        data = np.genfromtxt(
            input_star_file,
            skip_header=linenumber,
            dtype=dtype_list
            )

        return data, \
            np.array(header_names), np.array(input_star_file)

    else:
        return None

    assert(False)


def create_header_string(header_names):
    """Create and return the header string"""
    header_string = ['\ndata_\n\nloop_\n']
    for idx, entry in enumerate(header_names):
        header_string.append('{0} #{1}\n'.format(entry, idx+1))

    return ''.join(header_string)


def write_star_file(output_array, header_string, output_file, outliers, do_discard_outlier):
    """Write an output star file"""
    if do_discard_outlier:
        # Delete outliers
        not_outliers = np.where(outliers==0)[0]
    else:
        not_outliers = np.index_exp[:]

    # Write output
    with open('{0}'.format(output_file), 'w') as f:
        f.write(header_string)

        for row in output_array[not_outliers]:
            for element in row:
                if isinstance(element, (float, np.floating)):
                    text = '{:> 15.6f}'.format(element)
                elif isinstance(element, (int, np.integer)):
                    text = '{:> 7d}'.format(element)
                elif isinstance(element, (str, np.character)):
                    text = '{:>{}s}'.format(
                        element, len(element) + 6
                        )
                else:
                    assert False
                f.write(text)
            f.write('\n')
