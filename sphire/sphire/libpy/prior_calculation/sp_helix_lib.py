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
import os
import numpy.lib.recfunctions as rec
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import sp_helix_sphire as mhs
from . import sp_helix_relion as mhr
from . import sp_helix_prior as mhp


def identify_outliers(prior_tracker):
    """
    Identify outliers: if one angle is outlier, the whole filament is marked as an outlier.
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


def import_data_sphire(tracker, group_id, symclass, params_file=None, index_file=None):
    """
    Import the original stack information and create a tracker for the following calculations.
    :tracker: File name or dictionary - if dictionary it needs to contain the key stack or stack_prior:
    :params_file: SPHIRE - File containing the projection parameters:
    :index_file: SPHIRE - File containing the reference indices for the original stack:
    :return: prior_tracker containing the information from the sphire import:
    """

    prior_tracker = {}
    # Import the original stack for different input cases
    if isinstance(tracker, str):
        original_stack = mhs.import_sphire_stack(stack_path=tracker, group_id=group_id)
    elif isinstance(tracker, dict):
        # Only load the stack if it is not already loaded
        if 'stack_prior' in tracker['constants']:
            original_stack = tracker['constants']['stack_prior']
        else:
            original_stack = mhs.import_sphire_stack(stack_path=tracker['constants']['stack'], group_id=group_id)
    else:
        print('Unreachable code!')
        assert(False)

    # Import parameter and indices file and create a substack
    parameter = mhs.import_sphire_params(params_file, symclass=symclass)
    indices = mhs.import_sphire_index(index_file)
    substack = original_stack[indices['stack_idx']]

    # Combine the substack and the parameters information
    prior_tracker['array'] = rec.merge_arrays((substack, parameter, indices), usemask=False, flatten=True)

    # Names of the different array columns
    prior_tracker['order'] = 'source_n'

    prior_tracker['micrograph_id'] = 'ptcl_source_image'
    prior_tracker['data_n'] = 'data_n'

    prior_tracker[group_id] = group_id
    if group_id in original_stack.dtype.names and 'filament' not in group_id:
        prior_tracker['segment_id'] = "data_n"
    elif "filament" in original_stack.dtype.names and \
            "data_n" in original_stack.dtype.names and \
            'filament' in group_id:
        prior_tracker['segment_id'] = "data_n"
    elif "filament_id" in original_stack.dtype.names and \
            "segment_id" in original_stack.dtype.names and \
            "data_n" in original_stack.dtype.names and \
            'filament' in group_id:
        prior_tracker['segment_id'] = "segment_id"
    else:
        print ("ERROR: 'filament_id, segment_id and data_n' or 'filament and data_n' have to be present in the bdb header")
        print ("\t\t you should not be here because I prevent this kind of error in 'sp_helix_sphire.py' in the 'import_sphire_stack' function")
        exit(1)



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


def expand_and_order_array(prior_tracker, group_id):
    """
    Create an array with all columns present and sort it:
    micrograph_id > filament_id > segment_id
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
        prior_tracker[group_id],
        prior_tracker['segment_id']
        ]
    prior_tracker['array'] = np.sort(array_new, order=order)

    return prior_tracker


def loop_filaments(prior_tracker):
    """
    Loop over the filament to calculate prior values
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
