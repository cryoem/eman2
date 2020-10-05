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
import scipy.stats as so
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from . import sp_helix_lib as mhl

fit_index = 0

def wrapped_distribution(array):
    """Calculate a wrapped normal distribution"""
    nr_data = len(array)
    summe = np.sum(np.exp(np.radians(array)*1j))
    complex_mean = summe/float(nr_data)
    complex_angle = np.angle(complex_mean)

    R2 = np.real(complex_mean * np.conj(complex_mean))
    #R2_e = (nr_data/float(nr_data-1)) * (R2 - 1/float(nr_data))
    std = np.sqrt(np.log(1/float(R2)))

    mean_list = [np.degrees(complex_angle) for entry in range(nr_data)]
    std_list = [np.degrees(std) for entry in range(nr_data)]

    return mean_list, std_list


def get_filaments(prior_tracker, group_id):
    """Calculate the size and members of each filament with numpy"""
    # Create a list of unique filament names with indices
    filament_names, inverse = np.unique(
        prior_tracker['array'][prior_tracker[group_id]],
        return_inverse=True
        )
    # Loop through the indices and create subarrays
    filament_list = []
    for idx in range(len(filament_names)):
        indices_filament = np.where(inverse == idx)[0]
        filament_list.append(prior_tracker['array'][indices_filament])

    prior_tracker['array_filament'] = np.array(filament_list)

    return prior_tracker


def get_local_mean(input_array, rotate_angle, tolerance):
    """Calculate the local mean of the rotated array"""
    # Get the values of the outliers and non outliers
    inside_tolerance, outside_tolerance = get_tolerance_outliers(input_array, tolerance)

    # Calculate the local mean and rotate the coordinate system
    local_mean = np.mean(inside_tolerance)
    mean = subtract_and_adjust_angles(rotate_angle, -local_mean, 180, -180)
    subtract_and_adjust_angles(input_array, local_mean, 180, -180)

    return mean, len(outside_tolerance)


def rotate_angles_median(data_rotated, plot={'do_plot':False}):
    """Rotate angles to match 0 with the median"""
    if plot['do_plot']:
        mhl.plot_polar('adjusted_angles', data_rotated, 0, 180, -180, plot=plot)

    # Calculate where the most angles are
    diff_from_zero = np.median(np.abs(data_rotated))
    nr_positive = len(data_rotated[data_rotated >= 0])
    nr_negative = len(data_rotated) - nr_positive
    if nr_negative > nr_positive:
        current_median = -diff_from_zero
    else:
        current_median = diff_from_zero

    if plot['do_plot']:
        mhl.plot_polar('abs_values', np.abs(data_rotated), 0, 180, -180, plot=plot)
        mhl.plot_polar('abs_values_diff', np.abs(data_rotated), 0, 180, -180, mean=diff_from_zero, plot=plot)
        mhl.plot_polar('abs_values_median', data_rotated, 0, 180, -180, mean=current_median, plot=plot)

    iteration = 0
    rotate_angle = current_median
    current_median_old = current_median

    while iteration < 100:
        # In place subtraction
        subtract_and_adjust_angles(data_rotated, current_median, 180, -180)
        # Use the median for uneven numbers and the left median for even
        if len(data_rotated) % 2 == 0:
            median = np.sort(data_rotated)[len(data_rotated)//2]
        else:
            median = np.median(data_rotated)
        if plot['do_plot']:
            mhl.plot_polar('current_median_{0}'.format(iteration), data_rotated, rotate_angle, 180, -180, mean=median, old_mean=0, plot=plot)
        # Caluclate the new rotation angle
        rotate_angle = subtract_and_adjust_angles(rotate_angle, -median, 180, -180)

        # Check abort criteria
        if median == 0:
            break
        elif median == current_median_old:
            break
        else:
            current_median_old = current_median
            current_median = median
        iteration += 1

    if plot['do_plot']:
        mhl.plot_polar('rotated_data_median', data_rotated, rotate_angle, 180, -180, old_mean=0, plot=plot)

    return rotate_angle


def subtract_and_adjust_angles(data_rotated, value, angle_max, angle_min):
    """Subtract a value from an angle and adjust the angle range"""
    if isinstance(data_rotated, np.int64) or \
            isinstance(data_rotated, np.float64):
        data_subtracted = np.subtract(data_rotated, value)
        if data_subtracted > angle_max:
            data_subtracted -= 360
        elif data_subtracted <= angle_min:
            data_subtracted += 360
        else:
            pass
        return data_subtracted

    elif isinstance(data_rotated, np.ndarray):
        np.subtract(data_rotated, value, out=data_rotated)
        np.subtract(data_rotated, 360, out=data_rotated, where=data_rotated>angle_max)
        np.add(data_rotated, 360, out=data_rotated, where=data_rotated<angle_min)
        return None

    else:
        assert(False)


def rotate_angles_mean(data_rotated, rotate_angle, tol_mean, plot={'do_plot':False}):
    """Get filament outliers based on tolerance"""

    # Calculate the local mean
    rotate_angle, nr_outliers = get_local_mean(
        data_rotated,
        rotate_angle,
        tol_mean
        )
    rotate_angle_old = rotate_angle
    iteration = 0
    while iteration < 100:
        # Calculate the local mean
        rotate_angle, nr_outliers = get_local_mean(
            data_rotated,
            rotate_angle,
            tol_mean
            )

        # Abort criteria
        if rotate_angle == rotate_angle_old:
            break
        else:
            rotate_angle_old = rotate_angle
        iteration += 1

    return rotate_angle, nr_outliers


def identify_outliers_deg(data_rotated, tolerance, tolerance_filament, nr_outliers, plot={'do_plot':False}):
    inside_tolerance_idx, outside_tolerance_idx = find_tolerance_outliers(
        data_rotated, tolerance, nr_outliers
        )

    if len(outside_tolerance_idx) / float(len(data_rotated)) >= tolerance_filament:
        is_outlier = True
    else:
        is_outlier = False

    return is_outlier, inside_tolerance_idx, outside_tolerance_idx


def identify_outliers_std(std, tolerance, tolerance_std, plot={'do_plot':False}):
    if tolerance_std * std > tolerance:
        is_outlier = True
    else:
        is_outlier = False

    return is_outlier


def get_tolerance_outliers(input_array, tolerance):
    """Calculate the outliers"""
    inside_tolerance = input_array[np.logical_and(
        input_array >= -tolerance,
        input_array <= tolerance
        )]
    outside_tolerance = input_array[np.logical_or(
        input_array <= -tolerance,
        input_array >= tolerance
        )]

    return inside_tolerance, outside_tolerance


def find_tolerance_outliers(input_array, tolerance, nr_outliers):
    """Calculate the outliers"""
    if nr_outliers < len(input_array)/float(2):
        outside_tolerance = np.where(np.logical_or(
            input_array < -tolerance,
            input_array > tolerance
            ))[0]
        inside_tolerance = np.array(
            [i for i in range(len(input_array)) if i not in outside_tolerance]
            )
    else:
        inside_tolerance = np.where(np.logical_and(
            input_array >= -tolerance,
            input_array <= tolerance
            ))[0]
        outside_tolerance = np.array(
            [i for i in range(len(input_array)) if i not in inside_tolerance]
            )

    return inside_tolerance, outside_tolerance


def calculate_prior_values_running(data_rotated, prior_array, window_size, inside_tol_idx, outside_tol_idx, plot={'do_plot':False}):
    if window_size < 3:
        print("Window size needs to be at least 3: Changed window size to 3 and continue")
        window_size = 3

    # Check how to procede
    if len(inside_tol_idx) <= window_size * 1.5:
        # Fill the prior array with zeros
        prior_array.fill(0)
        return True
    else:
        # Calculate running averages
        assert(len(inside_tol_idx) > window_size * 1.5)

    # Calculate running averages
    mean_list = []
    mean_x_list = []
    number_of_means = len(data_rotated) - window_size + 1
    if window_size / float(2) % 2 == 0:
        min_x_mean = window_size // 2 - 1
    else:
        min_x_mean = window_size // 2
    max_x_mean = min_x_mean + number_of_means - 1

    inside_tol_idx = set(inside_tol_idx)
    for mean_idx in range(number_of_means):
        summation = 0
        skip = 0
        for number_idx in range(mean_idx, mean_idx + window_size):
            if number_idx in inside_tol_idx:
                summation += data_rotated[number_idx]
            else:
                assert(number_idx in outside_tol_idx)
                skip += 1
        if skip == window_size:
            continue
        else:
            mean_value = summation / float(window_size - skip)
            mean_list.append(mean_value)
            mean_x_list.append(min_x_mean + mean_idx)

    set_x_list = set(mean_x_list)
    for i in range(int(min_x_mean), int(max_x_mean)+1):
        if i not in set_x_list:
            idx = i - min_x_mean
        else:
            continue

        temp_mean = []
        if idx-1 < 0:
            pass
        else:
            temp_mean.append(mean_list[idx-1])
        if idx >= len(mean_x_list):
            pass
        else:
            temp_mean.append(mean_list[idx])
        mean_list.insert(idx, sum(temp_mean)/2)
        mean_x_list.insert(idx, i)

    min_value = mean_list[0]
    max_value = mean_list[-1]
    for i in reversed(range(min_x_mean)):
        mean_list.insert(0, min_value)
        mean_x_list.insert(0, i)
    for i in range(max_x_mean, len(data_rotated)-1):
        mean_list.append(max_value)
        mean_x_list.append(i)

    for idx in range(len(data_rotated)):
        prior_array[idx] = mean_list[idx]

    if plot['do_plot']:
        global fit_index
        plt.plot(mean_x_list, mean_list, 'x', label='running averages')
        plt.plot(range(len(data_rotated)), data_rotated, 'o', label='data')
        plt.ylim([-90, 90])
        plt.legend(loc='best')
        plt.xlabel('Particle helix id')
        plt.ylabel('Relative angle / degree')
        plt.grid()
        plt.savefig('{0}_fit_{1}.png'.format(plot['prefix'], fit_index))
        plt.clf()
        fit_index += 1

    return True


def calc_chi_square(y_values, x_values, params, fit_dim, func):
    points = len(y_values) - fit_dim - 1
    if points <= 1:
        points = 0.000000001
    else:
        pass
    return np.sum((y_values - func(x_values, *params))**2) / points


def calculate_prior_values_fit(data_rotated, prior_array, inside_tol_idx, plot={'do_plot':False}):
    func_linear = lambda x, a, b: a*x + b
    func_square = lambda x, a, b, c: a*x**2 + b*x + c
    func_cube = lambda x, a, b, c, d: a*x**3 + b*x**2 + c*x + d

    func_dict = {
        1: {'func': func_linear, 'params': None},
        2: {'func': func_square, 'params': None},
        3: {'func': func_cube, 'params': None},
        }

    chi_min = 99999999
    fit_min = 0
    for fit_dim in range(1, 4):
        try:
            params = np.polyfit(inside_tol_idx, data_rotated[inside_tol_idx], fit_dim)
        except ValueError:
            continue
        func = func_dict[fit_dim]['func']
        func_dict[fit_dim]['params'] = params

        chi_square = calc_chi_square(
            y_values=data_rotated[inside_tol_idx],
            x_values=inside_tol_idx,
            params=params,
            func=func,
            fit_dim=fit_dim
            )

        if chi_square < chi_min:
            chi_min = chi_square
            fit_min = fit_dim

    try:
        func_min = func_dict[fit_min]['func']
    except KeyError:
        for idx, entry in enumerate(data_rotated):
            prior_array[idx] = entry
    else:
        params_min = func_dict[fit_min]['params']
        for idx in range(len(data_rotated)):
            prior_array[idx] = func_min(idx, *params_min)

        if plot['do_plot']:
            global fit_index
            plt.plot(inside_tol_idx, data_rotated[inside_tol_idx], 'x', label='fitted data')
            x = np.linspace(0, len(data_rotated)-1, 10000)
            plt.plot(x, func_min(x, *params_min), label='regression')
            plt.legend(loc='best')
            plt.xlabel('Particle helix id')
            plt.ylabel('Relative angle / degree')
            plt.grid()
            plt.title('fit_dim: {0}'.format(fit_min))
            plt.savefig('{0}_fit_{1}.png'.format(plot['prefix'], fit_index))
            plt.clf()
            fit_index += 1


def mark_as_outlier(array, is_outlier, force_outlier, inside_tol_idx, outside_tol_idx):
    """
    Mark particles as outliers
    :array: Outlier array:
    :is_outlier: Variable, if the filament is an outlier:
    :force_outlier: Variable, if the particles are forced into the filament structure:
    :inside_tol_idx: Array of indices of non outliers:
    :outside_tol_idx: Array of indices of outliers:
    :return: None, do in place:
    """

    if is_outlier:
        array.fill(1)
    elif force_outlier:
        array.fill(0)
    else:
        array[inside_tol_idx] = 0
        array[outside_tol_idx] = 1


def calculate_prior_values_linear(data_rotated, prior_array, window_size, inside_tol_idx, outside_tol_idx, plot={'do_plot':False}):
    """Calculate running mean of the filament array"""

    if window_size < 3:
        print("Window size needs to be at least 3: Changed window size to 3 and continue")
        window_size = 3

    # Check how to procede
    if len(inside_tol_idx) <= window_size * 1.5:
        # Fill the prior array with zeros
        prior_array.fill(0)
        return True
    else:
        # Calculate running averages
        assert(len(inside_tol_idx) > window_size * 1.5)

    # Calculate running averages
    mean_list = []
    mean_x_list = []
    number_of_means = len(data_rotated) - window_size + 1
    if window_size / float(2) % 2 == 0:
        min_x_mean = window_size / float(2) - 0.5
    else:
        min_x_mean = window_size // float(2)

    inside_tol_idx = set(inside_tol_idx)
    for mean_idx in range(number_of_means):
        summation = 0
        skip = 0
        for number_idx in range(mean_idx, mean_idx + window_size):
            if number_idx in inside_tol_idx:
                summation += data_rotated[number_idx]
            else:
                assert(number_idx in outside_tol_idx)
                skip += 1
        if skip == window_size:
            continue
        else:
            mean_value = summation / float(window_size - skip)
            mean_list.append(mean_value)
            mean_x_list.append(min_x_mean + mean_idx)

    # Calculate linear regression
    mean_for_fit = np.array(mean_list)
    x_values_fit = np.array(mean_x_list)
    slope, intercept, rvalue, pvalue, stderr = so.linregress(x_values_fit, mean_for_fit)

    x_values_fit_mean = np.mean(x_values_fit)
    x_values_fit_variance = np.sum((x_values_fit - x_values_fit_mean)**2)
    slope_sd = stderr * np.sqrt(1/float(len(x_values_fit)) + x_values_fit_mean**2/x_values_fit_variance)
    intercept_sd = stderr * np.sqrt(1/x_values_fit_variance)

    for idx in range(len(data_rotated)):
        prior_array[idx] = slope * idx + intercept

    if plot['do_plot']:
        global fit_index
        plt.plot(x_values_fit, mean_for_fit, 'x', label='running averages')
        plt.plot(range(len(data_rotated)), data_rotated, 'o', label='data')
        x = np.linspace(0, len(data_rotated)-1, 10000)
        plt.plot(x, slope * x + intercept, label='linear regression')
        plt.ylim([-90, 90])
        plt.legend(loc='best')
        plt.xlabel('Particle helix id')
        plt.ylabel('Relative angle / degree')
        plt.grid()
        plt.title('stderr:{0}\nslopeerr:{1}\nintercepterr:{2}'.format(stderr, slope_sd, intercept_sd))
        plt.savefig('{0}_fit_{1}.png'.format(plot['prefix'], fit_index))
        plt.clf()
        fit_index += 1

    return True

