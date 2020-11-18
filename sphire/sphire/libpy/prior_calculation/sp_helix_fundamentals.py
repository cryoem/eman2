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

from . import sp_helix_prior as mhp
from . import sp_helix_lib as mhl


def calculate_priors(
        tracker,
        symclass,
        params_file=None,
        index_file=None,
        group_id=False, # Name of the column to check | E.g. filament_id, ISAC_class_id
        typ='sphire',
        tol_psi=30,
        tol_theta=15,
        tol_filament=0.2,
        tol_std=1,
        tol_mean=30,
        outlier_method='deg',
        prior_method='fit',
        force_outlier=False,
        plot=False,
        plot_lim=4,
        window_size=3,
        remove_outlier=False,
        node=0
        ):
    """
    Calculate prior values for the provided parameters

    prior_tracker keys:
    array
    order
    micrograph_id
    filament_id
    segment_id
    angle_names
    angle_max
    angle_min
    output_columns
    outlier
    output_file_params
    output_file_index
    output_dir
    idx_angle
    idx_angle_prior
    idx_angle_rot
    array_filament
    tol_filament
    """

    print('Prior calculation settings:')
    print('-> params_file:', params_file)
    print('-> index_file:', index_file)
    print('-> group_id:', group_id)
    print('-> typ:', typ)
    print('-> tol_psi:', tol_psi)
    print('-> tol_theta:', tol_theta)
    print('-> tol_filament:', tol_filament)
    print('-> tol_std:', tol_std)
    print('-> tol_mean:', tol_mean)
    print('-> outlier_method:', outlier_method)
    print('-> prior_method:', prior_method)
    print('-> force_outlier:', force_outlier)
    print('-> plot:', plot)
    print('-> plot_lim:', plot_lim)
    print('-> window_size:', window_size)
    print('-> remove_outlier:', remove_outlier)
    print('-> node:', node)

    # Import the stack and get parameters
    if typ == 'sphire':
        prior_tracker = mhl.import_data_sphire(
            tracker=tracker,
            index_file=index_file,
            params_file=params_file,
            group_id=group_id,
            symclass=symclass,
            )
    elif typ == 'relion':
        prior_tracker = mhl.import_data_relion(file_name=tracker)
    else:
        print('Unreachable code! Typ {0} not known! Supports sphire, relion'.format(
            typ
            ))
        return None

    # Create one huge array and sort the array
    prior_tracker = mhl.expand_and_order_array(prior_tracker=prior_tracker, group_id=group_id)

    # Seperate each filament into an array
    prior_tracker = mhp.get_filaments(prior_tracker=prior_tracker, group_id=group_id)

    # Add stuff to prior_tracker
    tolerance_list = [tol_theta, tol_psi]
    prior_tracker['plot'] = plot
    prior_tracker['plot_lim'] = plot_lim
    prior_tracker['window_size'] = window_size
    prior_tracker['tol_filament'] = tol_filament
    prior_tracker['tol_std'] = tol_std
    prior_tracker['tol_mean'] = tol_mean
    prior_tracker['node'] = node
    prior_tracker['apply_method'] = outlier_method
    prior_tracker['prior_method'] = prior_method
    prior_tracker['do_discard_outlier'] = remove_outlier
    prior_tracker['force_outlier'] = force_outlier
    # Execute calculation for each angle
    for idx, angle in enumerate(prior_tracker['angle_names']):
        prior_tracker['tolerance'] = tolerance_list[idx]
        prior_tracker['angle_prior'] = angle[prior_tracker['idx_angle_prior']]
        prior_tracker['angle_rot'] = angle[prior_tracker['idx_angle_rot']]
        prior_tracker['angle'] = angle[prior_tracker['idx_angle']]

        # Loop over the filament
        prior_tracker = mhl.loop_filaments(prior_tracker=prior_tracker)

    # Combine arrays and sort the combined array
    mhl.combine_and_order_filaments(prior_tracker=prior_tracker)

    # Print outliers
    IDX_ANGLE = prior_tracker['idx_angle']
    for angle in prior_tracker['angle_names']:
        column = 'outlier_{0}'.format(angle[IDX_ANGLE])
        print('==>', column, len(prior_tracker['array'][column][prior_tracker['array'][column] == 1]))

    # Identify outliers
    mhl.identify_outliers(prior_tracker=prior_tracker)

    # Write output
    if typ == 'sphire':
        mhl.export_data_sphire(prior_tracker=prior_tracker)
    elif typ == 'relion':
        mhl.export_data_relion(prior_tracker=prior_tracker)
    else:
        print('Unreachable code! Typ {0} not known! Supports sphire, relion'.format(
            typ
            ))
        return None

    return_list = [
        prior_tracker['array'][prior_tracker['outlier']],
        '{0}_not_applied.txt'.format('{0}_prior.txt'.format(prior_tracker['output_file_params'])),
        '{0}_not_applied.txt'.format('{0}_prior.txt'.format(prior_tracker['output_file_index'])),
        ]

    return return_list


if __name__ == '__main__':
    import shutil
    import os
    shutil.copy2('../tests/params_raw.txt', 'params.txt')
    shutil.copy2('../tests/index_raw.txt', 'index.txt')
    calculate_priors(
        'bdb:../tests/stack',
        params_file='params.txt',
        index_file='index.txt',
        typ='sphire',
        tol_psi=30,
        tol_theta=15,
        tol_filament=0.2,
        tol_std=1,
        tol_mean=30,
        outlier_method='deg',
        prior_method='fit',
        plot=False,
        plot_lim=4,
        window_size=3,
        remove_outlier=False,
        node=0
        )
    os.remove('params.txt')
    os.remove('params_prior.txt_not_applied.txt')
    os.remove('index.txt')
    os.remove('index_prior.txt_not_applied.txt')

