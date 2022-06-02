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
import sys
import EMAN2_cppwrap
import numpy as np
import shutil


def import_sphire_stack(stack_path, group_id):
    """Import the necessary data from a bdb/hdf stack"""
    dtype_list = [
        ('ptcl_source_image', '|U200'),
        ('filament', '|U200'),
        ('filament_id', '|U200'),
        ('segment_id', '<i8'),
        ('data_n', '<i8'),
        ('ISAC_class_id', '<i8')
        ]
    is_in_dtype = False
    for entry in dtype_list:
        if group_id == entry[0]:
            is_in_dtype = True
            break
    if not is_in_dtype:
        try:
            data = EMAN2_cppwrap.EMUtil.get_all_attributes(stack_path, group_id)
        except KeyError:
            print('Group_id', group_id, 'needs to be present in the stack header!')
            sys.exit(1)
        else:
            dtype = type(data)
        dtype_list.append((group_id, '|U200'))

    imported_data = []
    bad_idx = []
    for idx, entry in enumerate(dtype_list):
        try:
            data = EMAN2_cppwrap.EMUtil.get_all_attributes(stack_path, entry[0])
        except KeyError:
            bad_idx.append(idx)
            if entry[0] == group_id:
                print('Group_id', group_id, 'needs to be present in the stack header!')
                sys.exit(1)
        else:
            imported_data.append(data)
    for idx in reversed(bad_idx):
        del dtype_list[idx]

    data_array = np.empty(len(imported_data[0]), dtype=dtype_list)

    for idx, dtype in enumerate(dtype_list):
        data_array[dtype[0]] = imported_data[idx]

    return data_array


def import_sphire_params(input_file, symclass):
    """Import the params and index file"""

    dtype_import = [('phi', '<f8'), ('theta', '<f8'), ('psi', '<f8'), ('shift_x', '<f8'), ('shift_y', '<f8'), ('err1', '<f8'), ('err2', '<f8'), ('norm', '<f8')]
    dtype = dtype_import + [('source_n', '<i8')]
    #dtype = dtype_import + [('source_n', '<i8'), ('phi_old', '<f8'), ('theta_old', '<f8'), ('psi_old', '<f8')]

    data_import = np.genfromtxt(input_file, dtype=dtype_import)
    #reduced_angles = symclass.reduce_anglesets(data_import[['phi', 'theta', 'psi']].view(np.float64).reshape(data_import.shape + (-1,)), inc_mirror=0, tolistconv=False)

    #data['phi_old'] = data_import['phi']
    #data['theta_old'] = data_import['theta']
    #data['psi_old'] = data_import['psi']
    #data_import['phi'] = reduced_angles[:, 0]
    #data_import['theta'] = reduced_angles[:, 1]
    #data_import['psi'] = reduced_angles[:, 2]

    data = np.empty(len(data_import), dtype=dtype)
    data['source_n'] = np.arange(len(data))
    for name in data_import.dtype.names:
        data[name] = data_import[name]

    return data


def import_sphire_index(input_file):
    """Import the params and index file"""

    dtype = [('stack_idx', '<i8')]
    data = np.genfromtxt(input_file, dtype=dtype)

    return data


def write_params_file(array, names, file_name, file_name_old, prior_tracker):
    """Write sphire parameter file"""
    new_name_order = [
        'phi',
        'shift_x',
        'shift_y',
        'err1',
        'err2',
        'norm'
        ]
    for angle in prior_tracker['angle_names'][::-1]:
        new_name_order.insert(1, angle[prior_tracker['idx_angle_prior']])

    output_name = prepare_output(
        tracker=prior_tracker['tracker'],
        file_name=file_name,
        file_name_old=file_name_old
        )
    write_file(output_name=output_name, array=array, name_list=new_name_order, outlier_apply=prior_tracker['do_discard_outlier'], outlier_name=prior_tracker['outlier'])


def write_index_file(array, file_name, file_name_old, prior_tracker):
    """Write sphire index file"""

    output_name = prepare_output(
        tracker=prior_tracker['tracker'],
        file_name=file_name,
        file_name_old=file_name_old
        )
    write_file(output_name=output_name, array=array, name_list=['stack_idx'], outlier_apply=prior_tracker['do_discard_outlier'], outlier_name=prior_tracker['outlier'])

def write_file(output_name, array, name_list, outlier_apply, outlier_name):
    with open(output_name, 'w') as f:
        for element in array:
            if element[outlier_name] == 1:
                if outlier_apply:
                    continue
                else:
                    pass
            else:
                pass
            for name in name_list:
                if isinstance(element[name], (float, np.floating)):
                    text = '{:> 15.6f}'.format(element[name])
                elif isinstance(element[name], (int, np.integer)):
                    text = '{:> 7d}'.format(element[name])
                elif isinstance(element[name], (str, np.character)):
                    text = '{:>{}s}'.format(
                        element[name], len(element[name]) + 6
                        )
                else:
                    assert False
                f.write(text)
            f.write('\n')


def prepare_output(tracker, file_name, file_name_old):
    default_name = '{0}_not_applied.txt'.format(file_name)

    if isinstance(tracker, str):
        output_name = default_name
    elif isinstance(tracker, dict):
        #if tracker['constants']['apply_prior']:
        #    if tracker['state'] == 'RESTRICTED' or tracker['state'] == 'FINAL':
        #        shutil.move(file_name_old, file_name)
        #        output_name = file_name_old
        #    else:
        #        output_name = default_name

        #else:
        output_name = default_name
    else:
        print('Tracker instance "{0}" not known!'.format(type(tracker)))
        assert(False)

    return output_name
