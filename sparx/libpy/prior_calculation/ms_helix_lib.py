"""
author: markus stabrin 2017/07/26 (markus.stabrin@mpi-dortmund.mpg.de)

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
import os
import numpy.lib.recfunctions as rec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ms_helix_sphire as mhs
import ms_helix_relion as mhr
import ms_helix_prior as mhp


def identify_outliers(prior_tracker):
    """
    Identify outliers: if one angle is outlier, the whole filament is marked as an outlier.

    # No outliers
    >>> data = np.array([(0, 0, 0), (0, 0, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [0 0]

    # Filament 1 and Filament 2 outliers | Angle 1 and Angle 2 outlier
    >>> data = np.array([(1, 1, 0), (1, 1, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [1 1]

    # Filament 2 outliers | Angle 1 and Angle 2 outlier
    >>> data = np.array([(0, 0, 0), (1, 1, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [0 1]

    # Filament 1 outliers | Angle 1 and Angle 2 outlier
    >>> data = np.array([(1, 1, 0), (0, 0, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [1 0]

    # Filament 2 outliers | Angle 1 outlier
    >>> data = np.array([(0, 0, 0), (1, 0, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [0 1]

    # Filament 2 outliers | Angle 2 outlier
    >>> data = np.array([(0, 0, 0), (0, 1, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [0 1]

    # Filament 1 outliers | Angle 1 outlier
    >>> data = np.array([(1, 0, 0), (0, 0, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [1 0]

    # Filament 1 outliers | Angle 2 outlier
    >>> data = np.array([(0, 1, 0), (0, 0, 0)], dtype=[('outlier_angle1', '<i8'), ('outlier_angle2', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker = {'idx_angle': 0, 'array': data, 'angle_names': [['angle1'], ['angle2']]}
    >>> identify_outliers(prior_tracker)
    >>> print(prior_tracker['array']['outlier'])
    [1 0]

    :prior_tracker: dictionary containing the keys idx_angle, array, angle_names:
    :return: None, do in place calculation:
    """

    IDX_ANGLE = prior_tracker['idx_angle']
    for entry in prior_tracker['array']:
        is_outlier = False
        for angle in prior_tracker['angle_names']:
            column = 'outlier_{0}'.format(angle[IDX_ANGLE])
            if entry[column] == 1:
                is_outlier = True
            else:
                pass

        if is_outlier:
            entry['outlier'] = 1
        else:
            entry['outlier'] = 0


def combine_and_order_filaments(prior_tracker):
    """
    Combine filaments in one array and order the resulting array

    >>> data = np.array([(0, 0), (0, 0), (0, 0), (0, 0)], dtype=[('col1', '<i8'), ('col2', '<i8')])
    >>> data_filament1 = np.array([(1, 4), (2, 3)], dtype=[('col1', '<i8'), ('col2', '<i8')])
    >>> data_filament2 = np.array([(3, 2), (4, 1)], dtype=[('col1', '<i8'), ('col2', '<i8')])
    >>> data_filament = np.array([data_filament1, data_filament2])

    # Combine arrays | order col1 > col2
    >>> prior_tracker = {'array': data, 'array_filament': data_filament, 'order': ['col1', 'col2']}
    >>> combine_and_order_filaments(prior_tracker)
    >>> print(prior_tracker['array'][0])
    (1, 4)
    >>> print(prior_tracker['array'][1])
    (2, 3)
    >>> print(prior_tracker['array'][2])
    (3, 2)
    >>> print(prior_tracker['array'][3])
    (4, 1)

    # Combine arrays | order col2 > col1
    >>> prior_tracker = {'array': data, 'array_filament': data_filament, 'order': ['col2', 'col1']}
    >>> combine_and_order_filaments(prior_tracker)
    >>> print(prior_tracker['array'][0])
    (4, 1)
    >>> print(prior_tracker['array'][1])
    (3, 2)
    >>> print(prior_tracker['array'][2])
    (2, 3)
    >>> print(prior_tracker['array'][3])
    (1, 4)

    :prior_tracker: dictionary containing the keys array, array_filament, order:
    :return: prior_tracker containing the sorted combined array in the array slot:
    """
    array = np.empty(
        len(prior_tracker['array']), dtype=prior_tracker['array'].dtype.descr
        )
    index = 0
    for entry in prior_tracker['array_filament']:
        for row in entry:
            array[index] = row
            index += 1
    prior_tracker['array'] = np.sort(
        array, order=prior_tracker['order']
        )


def import_data_sphire(tracker, params_file=None, index_file=None):
    """
    Import the original stack information and create a tracker for the following calculations.

    >>> import shutil

    # Tracker is filename
    >>> tracker = 'bdb:../tests/stack'
    >>> params_file_raw = '../tests/index_raw.txt'
    >>> index_file_raw = '../tests/params_raw.txt'

    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> params_file = 'index.txt'
    >>> index_file = 'params.txt'
    >>> prior_tracker = import_data_sphire(tracker, params_file, index_file)
    >>> print(len(prior_tracker))
    16
    >>> print(len(prior_tracker['array']))
    684
    >>> print(prior_tracker.keys())
    ['angle_min', 'micrograph_id', 'segment_id', 'output_dir', 'output_file_params', 'idx_angle_prior', 'output_file_index', 'tracker', 'filament_id', 'angle_names', 'idx_angle', 'array', 'angle_max', 'order', 'idx_angle_rot', 'output_columns']
    >>> print(prior_tracker['array'][0])
    ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 0, 38.32375, 79.75458, 278.43793, 8.00109, 1.00109, 0.19693, 0.73653, 91259.15828, 0, 0)
    >>> print(prior_tracker['array'][-1])
    ('filt_Factin_ADP_cLys_0009_falcon2.hdf', 'filt_Factin_ADP_cLys_0009_falcon2.hdf0009', 7, 38.40719, 82.93683, 298.12543, -7.48667, 15.25721, 0.5188, 0.73356, 90890.99319, 683, 683)
    >>> print(prior_tracker['array'].dtype.names)
    ('ptcl_source_image', 'filament', 'data_n', 'phi', 'theta', 'psi', 'shift_x', 'shift_y', 'err1', 'err2', 'norm', 'source_n', 'stack_idx')
    >>> os.remove(params_file)
    >>> os.remove(index_file)

    # Tracker is dictionary and contains array
    >>> tracker = {'constants': {'stack_prior': prior_tracker['array'][['ptcl_source_image', 'filament', 'data_n']]}}
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> params_file = 'index.txt'
    >>> index_file = 'params.txt'
    >>> prior_tracker = import_data_sphire(tracker, params_file, index_file)
    >>> print(len(prior_tracker))
    16
    >>> print(len(prior_tracker['array']))
    684
    >>> print(prior_tracker.keys())
    ['angle_min', 'micrograph_id', 'segment_id', 'output_dir', 'output_file_params', 'idx_angle_prior', 'output_file_index', 'tracker', 'filament_id', 'angle_names', 'idx_angle', 'array', 'angle_max', 'order', 'idx_angle_rot', 'output_columns']
    >>> print(prior_tracker['array'][0])
    ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 0, 38.32375, 79.75458, 278.43793, 8.00109, 1.00109, 0.19693, 0.73653, 91259.15828, 0, 0)
    >>> print(prior_tracker['array'][-1])
    ('filt_Factin_ADP_cLys_0009_falcon2.hdf', 'filt_Factin_ADP_cLys_0009_falcon2.hdf0009', 7, 38.40719, 82.93683, 298.12543, -7.48667, 15.25721, 0.5188, 0.73356, 90890.99319, 683, 683)
    >>> print(prior_tracker['array'].dtype.names)
    ('ptcl_source_image', 'filament', 'data_n', 'phi', 'theta', 'psi', 'shift_x', 'shift_y', 'err1', 'err2', 'norm', 'source_n', 'stack_idx')
    >>> os.remove(params_file)
    >>> os.remove(index_file)

    # Tracker is dictionary and contains filename
    >>> tracker = {'constants': {'stack': 'bdb:../tests/stack'}}
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> params_file = 'index.txt'
    >>> index_file = 'params.txt'
    >>> prior_tracker = import_data_sphire(tracker, params_file, index_file)
    >>> print(len(prior_tracker))
    16
    >>> print(len(prior_tracker['array']))
    684
    >>> print(prior_tracker.keys())
    ['angle_min', 'micrograph_id', 'segment_id', 'output_dir', 'output_file_params', 'idx_angle_prior', 'output_file_index', 'tracker', 'filament_id', 'angle_names', 'idx_angle', 'array', 'angle_max', 'order', 'idx_angle_rot', 'output_columns']
    >>> print(prior_tracker['array'][0])
    ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 0, 38.32375, 79.75458, 278.43793, 8.00109, 1.00109, 0.19693, 0.73653, 91259.15828, 0, 0)
    >>> print(prior_tracker['array'][-1])
    ('filt_Factin_ADP_cLys_0009_falcon2.hdf', 'filt_Factin_ADP_cLys_0009_falcon2.hdf0009', 7, 38.40719, 82.93683, 298.12543, -7.48667, 15.25721, 0.5188, 0.73356, 90890.99319, 683, 683)
    >>> print(prior_tracker['array'].dtype.names)
    ('ptcl_source_image', 'filament', 'data_n', 'phi', 'theta', 'psi', 'shift_x', 'shift_y', 'err1', 'err2', 'norm', 'source_n', 'stack_idx')
    >>> os.remove(params_file)
    >>> os.remove(index_file)

    :tracker: File name or dictionary - if dictionary it needs to contain the key stack or stack_prior:
    :params_file: SPHIRE - File containing the projection parameters:
    :index_file: SPHIRE - File containing the reference indices for the original stack:
    :return: prior_tracker containing the information from the sphire import:
    """

    prior_tracker = {}
    # Import the original stack for different input cases
    if isinstance(tracker, basestring):
        original_stack = mhs.import_sphire_stack(stack_path=tracker)
    elif isinstance(tracker, dict):
        # Only load the stack if it is not already loaded
        if 'stack_prior' in tracker['constants']:
            original_stack = tracker['constants']['stack_prior']
        else:
            original_stack = mhs.import_sphire_stack(stack_path=tracker['constants']['stack'])
    else:
        print('Unreachable code!')
        assert(False)

    # Import parameter and indices file and create a substack
    parameter = mhs.import_sphire_params(params_file)
    indices = mhs.import_sphire_index(index_file)
    substack = original_stack[indices['stack_idx']]

    # Combine the substack and the parameters information
    prior_tracker['array'] = rec.merge_arrays((substack, parameter, indices), usemask=False, flatten=True)

    # Names of the different array columns
    prior_tracker['order'] = 'source_n'
    prior_tracker['micrograph_id'], prior_tracker['filament_id'], prior_tracker['segment_id'] = \
        original_stack.dtype.names
    prior_tracker['angle_names'] = [
        ['theta', 'theta_prior', 'theta_rot'],
        ['psi', 'psi_prior', 'psi_rot']
        ]
    prior_tracker['idx_angle'] = 0
    prior_tracker['idx_angle_prior'] = 1
    prior_tracker['idx_angle_rot'] = 2
    prior_tracker['angle_max'] = 360
    prior_tracker['angle_min'] = 0
    prior_tracker['output_columns'] = list(parameter.dtype.names)
    prior_tracker['output_file_params'] = params_file[:-len('.txt')]
    prior_tracker['output_file_index'] = index_file[:-len('.txt')]
    prior_tracker['output_dir'] = params_file[:-len(params_file.split('/')[-1])]
    if prior_tracker['output_dir'] == '':
        prior_tracker['output_dir'] = '.'

    prior_tracker['tracker'] = tracker

    return prior_tracker


def import_data_relion(file_name):
    """
    Import a relion star file

    # Tracker is filename
    >>> file_name = '../tests/data_test.star'
    >>> prior_tracker = import_data_relion(file_name=file_name)
    >>> print(len(prior_tracker))
    14
    >>> print(len(prior_tracker['array']))
    274
    >>> print(prior_tracker.keys())
    ['angle_min', 'micrograph_id', 'segment_id', 'output_dir', 'output_file', 'idx_angle_prior', 'filament_id', 'angle_names', 'idx_angle', 'array', 'angle_max', 'order', 'idx_angle_rot', 'output_columns']
    >>> print(prior_tracker['array'][0])
    (677.583313, 857.861694, 1, 0.0, '000001@Extract/job059/corrfull_1_39/Factin_ADP_cLys_0005_falcon2_DW.mrcs', 'corrfull_1_39/Factin_ADP_cLys_0005_falcon2_DW.mrc', 300.0, 20550.580078, 20305.132812, 23.187677, 0.0, 0.0, 1.0, 0.0, 0.1, 122807.0, 18.666667, 999.0, 0.006135, 1, -108.651069, 80.049875, -64.364309, -6.373347, 1.626653, 1, 0.808728, 1, 135071.7, 0.113652, 8, 0)
    >>> print(prior_tracker['array'][-1])
    (2896.461182, 1781.382568, 20, 592.000054, '000274@Extract/job059/corrfull_1_39/Factin_ADP_cLys_0005_falcon2_DW.mrcs', 'corrfull_1_39/Factin_ADP_cLys_0005_falcon2_DW.mrc', 300.0, 20550.580078, 20305.132812, 23.187677, 0.0, 0.0, 1.0, 0.0, 0.1, 122807.0, 18.666667, 999.0, 0.006135, 1, 149.725232, 87.340771, -138.469118, -4.901477, 2.098523, 1, 0.806169, 1, 134355.0, 0.172743, 20, 273)
    >>> print(prior_tracker['array'].dtype.names)
    ('_rlnCoordinateX', '_rlnCoordinateY', '_rlnHelicalTubeID', '_rlnHelicalTrackLength', '_rlnImageName', '_rlnMicrographName', '_rlnVoltage', '_rlnDefocusU', '_rlnDefocusV', '_rlnDefocusAngle', '_rlnSphericalAberration', '_rlnCtfBfactor', '_rlnCtfScalefactor', '_rlnPhaseShift', '_rlnAmplitudeContrast', '_rlnMagnification', '_rlnDetectorPixelSize', '_rlnCtfMaxResolution', '_rlnCtfFigureOfMerit', '_rlnGroupNumber', '_rlnAngleRot', '_rlnAngleTilt', '_rlnAnglePsi', '_rlnOriginX', '_rlnOriginY', '_rlnClassNumber', '_rlnNormCorrection', '_rlnRandomSubset', '_rlnLogLikeliContribution', '_rlnMaxValueProbDistribution', '_rlnNrOfSignificantSamples', 'source_n')

    :file_name: File name of the relion star file:
    :return: prior_tracker containing the information from the relion import:
    """
    prior_tracker = {}
    # Import the relion file
    array, header, path = mhr.import_star_file(input_star_file=file_name)

    # Names of the different array columns
    prior_tracker['order'] = 'source_n'
    prior_tracker['micrograph_id'] = '_rlnMicrographName'
    prior_tracker['filament_id'] = '_rlnHelicalTubeID'
    prior_tracker['segment_id'] = '_rlnImageName'
    prior_tracker['angle_names'] = [
        ['_rlnAngleTilt', '_rlnAngleTiltPrior', '_rlnAngleTiltRotated'],
        ['_rlnAnglePsi', '_rlnAnglePsiPrior', '_rlnAnglePsiRotated']
        ]
    prior_tracker['idx_angle'] = 0
    prior_tracker['idx_angle_prior'] = 1
    prior_tracker['idx_angle_rot'] = 2
    prior_tracker['angle_max'] = 180
    prior_tracker['angle_min'] = -180
    prior_tracker['output_columns'] = list(array.dtype.names)

    # Add a number list to the array
    array_stack_id = np.empty(len(array), dtype=[(prior_tracker['order'], '<i8')])
    array_stack_id[prior_tracker['order']] = np.arange(len(array))
    prior_tracker['array'] = rec.merge_arrays(
        (array, array_stack_id), usemask=False, flatten=True
        )
    prior_tracker['output_file'] = file_name[:-len('.star')]
    prior_tracker['output_dir'] = file_name[:-len(file_name.split('/')[-1])]
    if prior_tracker['output_dir'] == '':
        prior_tracker['output_dir'] = '.'

    return prior_tracker


def expand_and_order_array(prior_tracker):
    """
    Create an array with all columns present and sort it:
    micrograph_id > filament_id > segment_id

    >>> data = np.array([(0, 1, 2, 3, 4), (0, 1, 2, 3, 4), (0, 1, 2, 3, 4), (0, 1, 2, 3, 4)], dtype=[('col1', '<i8'), ('col2', '<i8'), ('col3', '<i8'), ('angle1', '<i8'), ('angle2', '<i8')])
    >>> prior_tracker = {}
    >>> prior_tracker['angle_names'] = [['angle1', 'angle1_prior', 'angle1_rot'], ['angle2', 'angle2_prior', 'angle2_rot']]
    >>> prior_tracker['array'] = data
    >>> prior_tracker['idx_angle'] = 0
    >>> prior_tracker['idx_angle_prior'] = 1
    >>> prior_tracker['idx_angle_rot'] = 2
    >>> prior_tracker['micrograph_id'] = 'col1'
    >>> prior_tracker['filament_id'] = 'col2'
    >>> prior_tracker['segment_id'] = 'col3'
    >>> prior_tracker['output_columns'] = ['col1', 'col2', 'col3']
    >>> prior_tracker = expand_and_order_array(prior_tracker=prior_tracker)
    >>> print(len(prior_tracker['output_columns']))
    5
    >>> print(len(prior_tracker['array']))
    4
    >>> print(len(prior_tracker['array'][0]))
    12

    :prior_tracker: Dictionary containing the keys idx_angle, idx_angle_prior, idx_angle_rot,
        array, angle_names, output_columns, micrograph_id, filament_id and segment_id:
    :return: prior_tracker containing the sorted combined array in the array slot:
    """

    # Add new angle names to the output dtype
    IDX_ANGLE = prior_tracker['idx_angle']
    IDX_PRIOR = prior_tracker['idx_angle_prior']
    IDX_ROT = prior_tracker['idx_angle_rot']
    dtype_new = prior_tracker['array'].dtype.descr
    for angle in prior_tracker['angle_names']:
        dtype_new.append((angle[IDX_PRIOR], '<f8'))
        dtype_new.append((angle[IDX_ROT], '<f8'))
        dtype_new.append(('outlier_{0}'.format(angle[IDX_ANGLE]), '<i8'))
        prior_tracker['output_columns'].append(angle[IDX_PRIOR])

    # Add outlier dtype
    prior_tracker['outlier'] = 'outlier'
    dtype_new.append((prior_tracker['outlier'], '<i8'))

    # Create a new large array that combines everything
    array_new = np.empty(len(prior_tracker['array']), dtype=dtype_new)
    for name in prior_tracker['array'].dtype.names:
        array_new[name] = prior_tracker['array'][name]
    for angle in prior_tracker['angle_names']:
        array_new[angle[IDX_ROT]] = np.copy(prior_tracker['array'][angle[IDX_ANGLE]])

    order = [
        prior_tracker['micrograph_id'],
        prior_tracker['filament_id'],
        prior_tracker['segment_id']
        ]
    prior_tracker['array'] = np.sort(array_new, order=order)

    return prior_tracker


def loop_filaments(prior_tracker):
    """
    Loop over the filament to calculate prior values

    >>> prior_tracker = {}
    >>> prior_tracker['array_filament'] = np.array([[ ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 0, 38.32375, 79.75458, 278.43793, 8.00109, 1.00109, 0.19693, 0.73653, 91259.15828, 0, 0, 6.042718322764362e-154, 79.75458, 2314885531238936880, 9.162199477420751e-72, 278.43793, 4051047449742946336, 3467807035425300512), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 1, 35.83593, 55.19391, 350.15674, -1.99891, -3.99891, 0.70511, 0.73476, 91040.31356, 1, 1, 4.948400661721011e+173, 55.19391, 8935143215481254912, 0.0, 350.15674, 8935145457454311936, 3487650908667907), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 2, 39.17223, 77.94179, 48.28169, -7.99891, -2.99891, 0.45521, 0.70617, 87497.71089, 2, 2, 0.0, 77.94179, 5157209670381677, 5.075883983234376e-116, 48.28169, 8317304086824383232, 8313495831334709343), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 3, 33.44401, 80.17834, 187.50038, -2.99891, -7.74279, 0.23672, 0.70772, 87689.01591, 3, 3, 0.0, 80.17834, 7595444411327411305, 7.698438237169649e+218, 187.50038, 1824520799039935077, 113044881408), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 4, 38.73295, 78.84809, 255.00044, 2.00109, 18.00109, 0.21688, 0.74511, 92321.562, 4, 4, 9.076859094468509e+223, 78.84809, 0, 0.0, 255.00044, 0, 0), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 5, 23.93829, 71.0277, 94.21904, -6.99891, 10.00109, 0.48858, 0.72281, 89558.54532, 5, 5, 0.0, 71.0277, 0, 0.0, 94.21904, 7161082258284898662, 6868064479304640884), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 6, 37.77144, 77.47843, 97.96918, -7.99891, -4.99891, 0.4357, 0.70789, 87710.85716, 6, 6, 0.0, 77.47843, 0, 0.0, 97.96918, 0, 0)]], dtype=[('ptcl_source_image', 'S200'), ('filament', 'S200'), ('data_n', '<i8'), ('phi', '<f8'), ('theta', '<f8'), ('psi', '<f8'), ('shift_x', '<f8'), ('shift_y', '<f8'), ('err1', '<f8'), ('err2', '<f8'), ('norm', '<f8'), ('source_n', '<i8'), ('stack_idx', '<i8'), ('theta_prior', '<f8'), ('theta_rot', '<f8'), ('outlier_theta', '<i8'), ('psi_prior', '<f8'), ('psi_rot', '<f8'), ('outlier_psi', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker['array'] = np.array([ ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 0, 38.32375, 79.75458, 278.43793, 8.00109, 1.00109, 0.19693, 0.73653, 91259.15828, 0, 0, 6.042718322764362e-154, 79.75458, 2314885531238936880, 9.162199477420751e-72, 278.43793, 4051047449742946336, 3467807035425300512), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 1, 35.83593, 55.19391, 350.15674, -1.99891, -3.99891, 0.70511, 0.73476, 91040.31356, 1, 1, 4.948400661721011e+173, 55.19391, 8935143215481254912, 0.0, 350.15674, 8935145457454311936, 3487650908667907), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 2, 39.17223, 77.94179, 48.28169, -7.99891, -2.99891, 0.45521, 0.70617, 87497.71089, 2, 2, 0.0, 77.94179, 5157209670381677, 5.075883983234376e-116, 48.28169, 8317304086824383232, 8313495831334709343), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 3, 33.44401, 80.17834, 187.50038, -2.99891, -7.74279, 0.23672, 0.70772, 87689.01591, 3, 3, 0.0, 80.17834, 7595444411327411305, 7.698438237169649e+218, 187.50038, 1824520799039935077, 113044881408), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 4, 38.73295, 78.84809, 255.00044, 2.00109, 18.00109, 0.21688, 0.74511, 92321.562, 4, 4, 9.076859094468509e+223, 78.84809, 0, 0.0, 255.00044, 0, 0), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 5, 23.93829, 71.0277, 94.21904, -6.99891, 10.00109, 0.48858, 0.72281, 89558.54532, 5, 5, 0.0, 71.0277, 0, 0.0, 94.21904, 7161082258284898662, 6868064479304640884), ('filt_Factin_ADP_cLys_0005_falcon2.hdf', 'filt_Factin_ADP_cLys_0005_falcon2.hdf0000', 6, 37.77144, 77.47843, 97.96918, -7.99891, -4.99891, 0.4357, 0.70789, 87710.85716, 6, 6, 0.0, 77.47843, 0, 0.0, 97.96918, 0, 0)], dtype=[('ptcl_source_image', 'S200'), ('filament', 'S200'), ('data_n', '<i8'), ('phi', '<f8'), ('theta', '<f8'), ('psi', '<f8'), ('shift_x', '<f8'), ('shift_y', '<f8'), ('err1', '<f8'), ('err2', '<f8'), ('norm', '<f8'), ('source_n', '<i8'), ('stack_idx', '<i8'), ('theta_prior', '<f8'), ('theta_rot', '<f8'), ('outlier_theta', '<i8'), ('psi_prior', '<f8'), ('psi_rot', '<f8'), ('outlier_psi', '<i8'), ('outlier', '<i8')])
    >>> prior_tracker['do_discard_outlier'] = False
    >>> prior_tracker['plot'] = False
    >>> prior_tracker['window_size'] = 3
    >>> prior_tracker['angle'] = 'theta'
    >>> prior_tracker['prior_method'] = 'fit'
    >>> prior_tracker['output_file_params'] = 'params'
    >>> prior_tracker['apply_method'] = 'deg'
    >>> prior_tracker['plot_lim'] = 4
    >>> prior_tracker['tracker'] = 'bdb:../tests/stack_small'
    >>> prior_tracker['output_dir'] = '.'
    >>> prior_tracker['angle_names'] = [['theta', 'theta_prior', 'theta_rot'], ['psi', 'psi_prior', 'psi_rot']]
    >>> prior_tracker['idx_angle'] = 0
    >>> prior_tracker['tol_filament'] = 0.2
    >>> prior_tracker['angle_prior'] = 'theta_prior'
    >>> prior_tracker['tolerance'] = 15
    >>> prior_tracker['node'] = 0
    >>> prior_tracker['angle_min'] = 0
    >>> prior_tracker['outlier'] = 'outlier'
    >>> prior_tracker['angle_rot'] = 'theta_rot'
    >>> prior_tracker['idx_angle_prior'] = 1
    >>> prior_tracker['output_file_index'] = 'index'
    >>> prior_tracker['angle_max'] = 360
    >>> prior_tracker['output_columns'] = ['phi', 'theta', 'psi', 'shift_x', 'shift_y', 'err1', 'err2', 'norm', 'source_n', 'theta_prior', 'psi_prior']
    >>> prior_tracker['tol_std'] = 1
    >>> prior_tracker['micrograph_id'] = 'ptcl_source_image'
    >>> prior_tracker['segment_id'] = 'data_n'
    >>> prior_tracker['tol_mean'] = 30
    >>> prior_tracker['filament_id'] = 'filament'
    >>> prior_tracker['order'] = 'source_n'
    >>> prior_tracker['idx_angle_rot'] = 2
    >>> prior_tracker['force_outlier'] = True
    >>> prior_tracker = loop_filaments(prior_tracker=prior_tracker)
    >>> print(isinstance(prior_tracker, dict))
    True
    """
    angle_prior = prior_tracker['angle_prior']
    angle_rot = prior_tracker['angle_rot']
    angle = prior_tracker['angle']
    window_size = prior_tracker['window_size']
    plot_dict = {}
    plot_dict['output_dir'] = '{0}/prior_images_{1}'.format(
        prior_tracker['output_dir'], prior_tracker['node']
        )
    if not os.path.exists(plot_dict['output_dir']):
        os.mkdir(plot_dict['output_dir'])

    for idx, filament in enumerate(prior_tracker['array_filament']):
        # Plot settings
        if idx < prior_tracker['plot_lim'] and prior_tracker['plot']:
            plot_dict['do_plot'] = True
            plot_dict['prefix'] = '{0}/{1}_{2}'.format(plot_dict['output_dir'], angle, idx)
        else:
            plot_dict['do_plot'] = False

        # Shift the angle_range to [180,-180]
        mhp.subtract_and_adjust_angles(
            data_rotated=filament[angle_rot],
            value=0,
            angle_max=180,
            angle_min=-180
            )

        # Rotate the angle range, so that the median is the new center
        rotate_angle = mhp.rotate_angles_median(
            data_rotated=filament[angle_rot],
            plot=plot_dict
            )

        # Rotate the angle range, so that the mean is the new center
        rotate_angle, nr_outliers = mhp.rotate_angles_mean(
            data_rotated=filament[angle_rot],
            rotate_angle=rotate_angle,
            tol_mean=prior_tracker['tol_mean'],
            plot=plot_dict
            )

        if plot_dict['do_plot']:
            plot_polar('rotated_data_mean', filament[angle_rot], rotate_angle, prior_tracker['angle_max'], prior_tracker['angle_min'], plot=plot_dict)

        # Calculate outliers based on the method
        if prior_tracker['apply_method'] == 'deg':
            # Identify based on degree
            is_outlier, inside_tol_idx, outside_tol_idx = mhp.identify_outliers_deg(
                data_rotated=filament[angle_rot],
                tolerance=prior_tracker['tolerance'],
                tolerance_filament=prior_tracker['tol_filament'],
                nr_outliers=nr_outliers,
                plot=plot_dict
                )
            if plot_dict['do_plot']:
                plot_polar('outliers_deg', filament[angle_rot], rotate_angle, prior_tracker['angle_max'], prior_tracker['angle_min'], plot=plot_dict, mean=0, tol=prior_tracker['tolerance'])

        elif prior_tracker['apply_method'] == 'std':
            is_outlier, inside_tol_idx, outside_tol_idx = mhp.identify_outliers_deg(
                data_rotated=filament[angle_rot],
                tolerance=90,
                tolerance_filament=prior_tracker['tol_filament'],
                nr_outliers=nr_outliers,
                plot=plot_dict
                )
            if plot_dict['do_plot']:
                plot_polar('outliers_deg', filament[angle_rot], rotate_angle, prior_tracker['angle_max'], prior_tracker['angle_min'], plot=plot_dict, mean=0, tol=90)

            if is_outlier:
                pass
            else:
                # Calculate standard deviation of angular distribution
                mean_list, std_list = mhp.wrapped_distribution(array=filament[angle][inside_tol_idx])
                # Calculate outliers based on std
                is_outlier = mhp.identify_outliers_std(
                    std=std_list[0],
                    tolerance=prior_tracker['tolerance'],
                    tolerance_std=prior_tracker['tol_std'],
                    plot=plot_dict
                    )
                if plot_dict['do_plot']:
                    plot_polar('outliers_std', filament[angle_rot], rotate_angle, prior_tracker['angle_max'], prior_tracker['angle_min'], plot=plot_dict, mean=0, tol=prior_tracker['tol_std']*std_list[0])
        else:
            print('\n', 'Apply method {0} not known! Use deg'.format(prior_tracker['apply_method']), '\n')
            is_outlier, inside_tol_idx, outside_tol_idx = mhp.identify_outliers_deg(
                data_rotated=filament[angle_rot],
                tolerance=prior_tracker['tolerance'],
                tolerance_filament=['tol_filament'],
                nr_outliers=nr_outliers,
                plot=plot_dict
                )

        # Mark as outlier
        mhp.mark_as_outlier(
            array=filament['outlier_{0}'.format(angle)],
            is_outlier=is_outlier,
            force_outlier=prior_tracker['force_outlier'],
            inside_tol_idx=inside_tol_idx,
            outside_tol_idx=outside_tol_idx
            )

        # Calculate prior values
        if is_outlier:
            filament[angle_prior] = filament[angle_rot]
        elif prior_tracker['prior_method'] == 'linear':
            mhp.calculate_prior_values_linear(
                data_rotated=filament[angle_rot],
                prior_array=filament[angle_prior],
                window_size=window_size,
                inside_tol_idx=inside_tol_idx,
                outside_tol_idx=outside_tol_idx,
                plot=plot_dict
                )
        elif prior_tracker['prior_method'] == 'fit':
            mhp.calculate_prior_values_fit(
                data_rotated=filament[angle_rot],
                prior_array=filament[angle_prior],
                inside_tol_idx=inside_tol_idx,
                plot=plot_dict
                )
        elif prior_tracker['prior_method'] == 'running':
            mhp.calculate_prior_values_running(
                data_rotated=filament[angle_rot],
                prior_array=filament[angle_prior],
                window_size=window_size,
                inside_tol_idx=inside_tol_idx,
                outside_tol_idx=outside_tol_idx,
                plot=plot_dict
                )
        else:
            print('prior_method', prior_tracker['prior_method'], 'not known! Use fit now')
            mhp.calculate_prior_values_fit(
                data_rotated=filament[angle_rot],
                prior_array=filament[angle_prior],
                inside_tol_idx=inside_tol_idx,
                plot=plot_dict
                )

        # Adjust angle range
        mhp.subtract_and_adjust_angles(
            filament[angle_prior],
            -rotate_angle,
            prior_tracker['angle_max'],
            prior_tracker['angle_min']
            )

        if plot_dict['do_plot']:
            plot_polar('final_angles', filament[angle_prior], 0, prior_tracker['angle_max'], prior_tracker['angle_min'], plot=plot_dict)

    return prior_tracker


def export_data_relion(prior_tracker):
    """
    Export the calculated priors for relion

    # Export without outliers | Remove outlier from file
    >>> data = np.array([(0, 1, 2, 3, 0), (5, 6, 7, 8, 0), (10, 11, 12, 13, 0), (15, 16, 17, 18, 0)], dtype=[('col1', '<i8'), ('col2', '<i8'), ('col3', '<i8'), ('angle1', '<i8'), ('outlier', '<i8')])
    >>> output_columns = ['col1', 'col2', 'col3']
    >>> output_file_name = 'doctest'
    >>> prior_tracker = {'array': data, 'output_columns': output_columns, 'output_file': output_file_name, 'outlier': 'outlier', 'do_discard_outlier': True}
    >>> export_data_relion(prior_tracker=prior_tracker)
    >>> output_file = 'doctest_prior.star'
    >>> os.path.exists(output_file)
    True
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    11
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(lines[-1].rstrip('\\n'))
         15     16     17
    >>> os.remove(output_file) 

    # Test with one outlier in the beginning | Remove outlier from file
    >>> data = np.array([(0, 1, 2, 3, 1), (5, 6, 7, 8, 0), (10, 11, 12, 13, 0), (15, 16, 17, 18, 0)], dtype=[('col1', '<i8'), ('col2', '<i8'), ('col3', '<i8'), ('angle1', '<i8'), ('outlier', '<i8')])
    >>> output_columns = ['col1', 'col2', 'col3']
    >>> output_file_name = 'doctest'
    >>> prior_tracker = {'array': data, 'output_columns': output_columns, 'output_file': output_file_name, 'outlier': 'outlier', 'do_discard_outlier': True}
    >>> export_data_relion(prior_tracker=prior_tracker)
    >>> output_file = 'doctest_prior.star'
    >>> os.path.exists(output_file)
    True
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    10
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(lines[-1].rstrip('\\n'))
         15     16     17
    >>> os.remove(output_file) 

    # Test with one outlier in the beginning and one in the end | Remove outlier from file
    >>> data = np.array([(0, 1, 2, 3, 1), (5, 6, 7, 8, 0), (10, 11, 12, 13, 0), (15, 16, 17, 18, 1)], dtype=[('col1', '<i8'), ('col2', '<i8'), ('col3', '<i8'), ('angle1', '<i8'), ('outlier', '<i8')])
    >>> output_columns = ['col1', 'col2', 'col3']
    >>> output_file_name = 'doctest'
    >>> prior_tracker = {'array': data, 'output_columns': output_columns, 'output_file': output_file_name, 'outlier': 'outlier', 'do_discard_outlier': True}
    >>> export_data_relion(prior_tracker=prior_tracker)
    >>> output_file = 'doctest_prior.star'
    >>> os.path.exists(output_file)
    True
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    9
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(lines[-1].rstrip('\\n'))
         10     11     12
    >>> os.remove(output_file) 

    # Test with one outlier in the beginning and one in the end | Dont remove outliers from file
    >>> data = np.array(
    ...     [(0, 1, 2, 3, 1), (5, 6, 7, 8, 0), (10, 11, 12, 13, 0), (15, 16, 17, 18, 1)],
    ...     dtype=[('col1', '<i8'), ('col2', '<i8'), ('col3', '<i8'), ('angle1', '<i8'), ('outlier', '<i8')]
    ...     )
    >>> output_columns = ['col1', 'col2', 'col3']
    >>> output_file_name = 'doctest'
    >>> prior_tracker = {'array': data, 'output_columns': output_columns, 'output_file': output_file_name, 'outlier': 'outlier', 'do_discard_outlier': False}
    >>> export_data_relion(prior_tracker=prior_tracker)
    >>> output_file = 'doctest_prior.star'
    >>> os.path.exists(output_file)
    True
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    11
    >>> with open(output_file, 'r') as read:
    ...     lines = read.readlines()
    ...     print(lines[-1].rstrip('\\n'))
         15     16     17
    >>> os.remove(output_file) 

    :prior_tracker: Dictionary containing the keys array, output_columns, output_file:
    :return: None, just write files
    """

    header_string = mhr.create_header_string(prior_tracker['output_columns'])
    mhr.write_star_file(
        output_array=prior_tracker['array'][prior_tracker['output_columns']],
        header_string=header_string,
        output_file='{0}_prior.star'.format(prior_tracker['output_file']),
        outliers=prior_tracker['array'][prior_tracker['outlier']],
        do_discard_outlier=prior_tracker['do_discard_outlier']
        )


def export_data_sphire(prior_tracker):
    """
    Export the calculated priors for sphire
    >>> import shutil

    >>> params_file_raw = '../tests/params_raw.txt'
    >>> index_file_raw = '../tests/index_raw.txt'
    >>> data = np.array(
    ...     [(0, 1, 2, 3, 0, 1, 2, 0, 0), (5, 6, 7, 8, 0, 1, 2, 0, 1), (10, 11, 12, 13, 0, 1, 2, 0, 2), (15, 16, 17, 18, 0, 1, 2, 1, 3)],
    ...     dtype=[('angle1_prior', '<i8'), ('phi', '<i8'), ('shift_x', '<i8'), ('shift_y', '<i8'), ('err1', '<i8'), ('err2', '<i8'), ('norm', '<i8'), ('outlier', '<i8'), ('stack_idx', '<i8')]
    ...     )
    >>> output_columns = ['col1', 'col2', 'col3']
    >>> output_file_params = 'params'
    >>> output_file_index = 'index'

    # One outlier | tracker dict RESTRICTED | Dont discard outliers
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> prior_tracker = {
    ...     'array': data,
    ...     'output_columns': output_columns,
    ...     'output_file_params': output_file_params,
    ...     'output_file_index': output_file_index,
    ...     'idx_angle_prior':1,
    ...     'angle_names': [['angle1', 'angle1_prior', 'angle1_rot']],
    ...     'outlier': 'outlier',
    ...     'tracker': {'constants': {'apply_prior': True}, 'state': 'RESTRICTED'},
    ...     'do_discard_outlier': False
    ...     }
    >>> export_data_sphire(prior_tracker=prior_tracker)
    >>> with open('{0}.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    4
         16     15     17     18      0      1      2
    >>> with open('{0}_prior.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
          38.40719      82.93683     298.12543      -7.48667      15.25721       0.51880       0.73356   90890.99319
    >>> with open('{0}.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    4
          3
    >>> with open('{0}_prior.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
    683
    >>> os.remove('{0}.txt'.format(output_file_index))
    >>> os.remove('{0}_prior.txt'.format(output_file_index))
    >>> os.remove('{0}.txt'.format(output_file_params))
    >>> os.remove('{0}_prior.txt'.format(output_file_params))

    # One outlier | tracker dict outlier RESTRICTED | Discard outliers
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> prior_tracker = {
    ...     'array': data,
    ...     'output_columns': output_columns,
    ...     'output_file_params': output_file_params,
    ...     'output_file_index': output_file_index,
    ...     'idx_angle_prior':1,
    ...     'angle_names': [['angle1', 'angle1_prior', 'angle1_rot']],
    ...     'outlier': 'outlier',
    ...     'tracker': {'constants': {'apply_prior': True}, 'state': 'RESTRICTED'},
    ...     'do_discard_outlier': True
    ...     }
    >>> export_data_sphire(prior_tracker=prior_tracker)
    >>> with open('{0}.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
         11     10     12     13      0      1      2
    >>> with open('{0}_prior.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
          38.40719      82.93683     298.12543      -7.48667      15.25721       0.51880       0.73356   90890.99319
    >>> with open('{0}.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
          2
    >>> with open('{0}_prior.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
    683
    >>> os.remove('{0}.txt'.format(output_file_index))
    >>> os.remove('{0}_prior.txt'.format(output_file_index))
    >>> os.remove('{0}.txt'.format(output_file_params))
    >>> os.remove('{0}_prior.txt'.format(output_file_params))

    # One outlier | tracker dict outlier EXHAUSTIVE | Discard outliers
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> prior_tracker = {
    ...     'array': data,
    ...     'output_columns': output_columns,
    ...     'output_file_params': output_file_params,
    ...     'output_file_index': output_file_index,
    ...     'idx_angle_prior':1,
    ...     'angle_names': [['angle1', 'angle1_prior', 'angle1_rot']],
    ...     'outlier': 'outlier',
    ...     'tracker': {'constants': {'apply_prior': True}, 'state': 'EXHAUSTIVE'},
    ...     'do_discard_outlier': True
    ...     }
    >>> export_data_sphire(prior_tracker=prior_tracker)
    >>> with open('{0}.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
          38.40719      82.93683     298.12543      -7.48667      15.25721       0.51880       0.73356   90890.99319
    >>> with open('{0}_prior.txt_not_applied.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
         11     10     12     13      0      1      2
    >>> with open('{0}.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
    683
    >>> with open('{0}_prior.txt_not_applied.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
          2
    >>> os.remove('{0}.txt'.format(output_file_index))
    >>> os.remove('{0}_prior.txt_not_applied.txt'.format(output_file_index))
    >>> os.remove('{0}.txt'.format(output_file_params))
    >>> os.remove('{0}_prior.txt_not_applied.txt'.format(output_file_params))

    # One outliers | tracker dict outlier EXHAUSTIVE | Discard outliers
    >>> shutil.copy2(index_file_raw, 'index.txt')
    >>> shutil.copy2(params_file_raw, 'params.txt')
    >>> prior_tracker = {
    ...     'array': data,
    ...     'output_columns': output_columns,
    ...     'output_file_params': output_file_params,
    ...     'output_file_index': output_file_index,
    ...     'idx_angle_prior':1,
    ...     'angle_names': [['angle1', 'angle1_prior', 'angle1_rot']],
    ...     'outlier': 'outlier',
    ...     'tracker': {'constants': {'apply_prior': True}, 'state': 'EXHAUSTIVE'},
    ...     'do_discard_outlier': True
    ...     }
    >>> export_data_sphire(prior_tracker=prior_tracker)
    >>> with open('{0}.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
          38.40719      82.93683     298.12543      -7.48667      15.25721       0.51880       0.73356   90890.99319
    >>> with open('{0}_prior.txt_not_applied.txt'.format(output_file_params), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
         11     10     12     13      0      1      2
    >>> with open('{0}.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    684
    683
    >>> with open('{0}_prior.txt_not_applied.txt'.format(output_file_index), 'r') as read:
    ...     lines = read.readlines()
    ...     print(len(lines))
    ...     print(lines[-1].rstrip('\\n'))
    3
          2
    >>> os.remove('{0}.txt'.format(output_file_index))
    >>> os.remove('{0}_prior.txt_not_applied.txt'.format(output_file_index))
    >>> os.remove('{0}.txt'.format(output_file_params))
    >>> os.remove('{0}_prior.txt_not_applied.txt'.format(output_file_params))

    :prior_tracker: Dictionary containing the keys array, output_columns, output_file_params, output_file_index :
    :return: None, just write files
    """

    mhs.write_params_file(
        array=prior_tracker['array'],
        names=prior_tracker['output_columns'],
        file_name='{0}_prior.txt'.format(prior_tracker['output_file_params']),
        file_name_old='{0}.txt'.format(prior_tracker['output_file_params']),
        prior_tracker=prior_tracker
        )
    mhs.write_index_file(
        array=prior_tracker['array'],
        file_name='{0}_prior.txt'.format(prior_tracker['output_file_index']),
        file_name_old='{0}.txt'.format(prior_tracker['output_file_index']),
        prior_tracker=prior_tracker
        )


plot_index = 0
def plot_polar(name, array, angle_rotation, angle_max, angle_min, mean=None, tol=None, old_mean=None, plot={'prefix': 'DEFAULT'}):
    """
    Do a polar plot

    # Test default values | angle_max [0, 360)
    >>> plot_polar(
    ...     name='doctest',
    ...     array=np.array([0, 1, 2, 3]),
    ...     angle_rotation=3,
    ...     angle_max=360,
    ...     angle_min=0
    ...     )
    >>> os.path.exists('DEFAULT_0_doctest.png')
    True
    >>> os.remove('DEFAULT_0_doctest.png')

    # Test default values | angle (-180, 180]
    >>> plot_polar(
    ...     name='doctest',
    ...     array=np.array([0, 1, 2, 3]),
    ...     angle_rotation=3,
    ...     angle_max=180,
    ...     angle_min=-180
    ...     )
    >>> os.path.exists('DEFAULT_1_doctest.png')
    True
    >>> os.remove('DEFAULT_1_doctest.png')

    # Test default values | angle [0, 90]
    >>> plot_polar(
    ...     name='doctest',
    ...     array=np.array([0, 1, 2, 3]),
    ...     angle_rotation=3,
    ...     angle_max=90,
    ...     angle_min=0
    ...     )
    >>> os.path.exists('DEFAULT_2_doctest.png')
    True
    >>> os.remove('DEFAULT_2_doctest.png')

    # Test no default values | angle_max [0, 360)
    >>> plot_polar(
    ...     name='doctest',
    ...     array=np.array([0, 1, 2, 3]),
    ...     angle_rotation=3,
    ...     angle_max=360,
    ...     angle_min=0,
    ...     mean=20,
    ...     tol=5,
    ...     old_mean=8
    ...     )
    >>> os.path.exists('DEFAULT_3_doctest.png')
    True
    >>> os.remove('DEFAULT_3_doctest.png')

    :name: Name of the plot:
    :array: Array to plot:
    :angle_rotation: Specifies the zero position of the plot:
    :angle_max: Maximum angle lable:
    :angle_min: Minimum angle lable:
    :mean: Value of the mean:
    :tol: Value of the tolerance value:
    :old_mean: Value of the old mean:
    :plot: Plot dictionary; Need to have the key prefix:
    :return: None, writes image to disc:
    """
    global plot_index
    # radar green, solid grid lines
    plt.rc('grid', color='#316931', linewidth=1, linestyle='-')
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # force square figure and square axes looks better for polar, IMO
    width, height = matplotlib.rcParams['figure.figsize']
    size = min(width, height)

    # make a square figure
    fig = plt.figure(figsize=(size, size))
    try:
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, axisbg='#d5de9c')
    except AttributeError:
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True, facecolor='#d5de9c')

    # Change the labels
    ax.set_yticklabels([])
    ax.set_theta_direction(-1)
    offset = -np.radians(angle_rotation) + np.pi/2
    ax.set_theta_offset(offset)
    if angle_max == 360:
        labels = [str(item) if item != 360 else str(0) for item in np.linspace(0, 360, 9) if item >= angle_min]
    elif angle_max == 180:
        labels = [str(item-360) if item > 180 else str(item) for item in np.linspace(0, 360, 9) if item >= angle_min]
    else:
        labels = [str(item) if item != 360 else str(0) for item in np.linspace(angle_min, angle_max, 9) if item >= angle_min]
    ax.set_xticklabels(labels)

    # Plot the data
    color = plt.cm.Dark2(np.linspace(0, 1, len(array)))
    for idx in range(len(array)):
        plt.arrow(np.radians(array[idx]), 0.01, 0, 0.9, alpha=1, width=0.005, edgecolor=color[idx], facecolor=color[idx], lw=2, zorder=5)

    if mean is not None:
        plt.arrow(np.radians(mean), 0.01, 0, 0.9, alpha=1, width=0.005, edgecolor='black', facecolor='black', lw=4, zorder=5)

    if old_mean is not None:
        plt.arrow(np.radians(old_mean), 0.01, 0, 0.7, alpha=1, width=0.005, edgecolor='blue', facecolor='blue', lw=4, zorder=5)

    if tol is not None:
        plt.arrow(np.radians(mean+tol), 0.01, 0, 0.5, alpha=0.5, width=0.005, edgecolor='black', facecolor='black', lw=2, zorder=5)
        plt.arrow(np.radians(mean-tol), 0.01, 0, 0.5, alpha=0.5, width=0.005, edgecolor='black', facecolor='black', lw=2, zorder=5)

    # Beautify
    plt.title(name)
    plt.savefig('{0}_{1}_{2}.png'.format(plot['prefix'], plot_index, name))
    plot_index += 1
    plt.close(fig)
