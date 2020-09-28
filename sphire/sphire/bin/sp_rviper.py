#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#
#

import EMAN2_cppwrap
import copy
import json
import matplotlib
import matplotlib.pyplot
import mpi
import optparse
import os
import random
import shutil
import six
from ..libpy import sp_applications
from ..libpy import sp_fundamentals
from ..libpy import sp_global_def
from ..libpy import sp_logger
from ..libpy import sp_multi_shc
from ..libpy import sp_statistics
from ..libpy import sp_user_functions
from ..libpy import sp_utilities
import string
import sys
import time
import itertools
import re
from future import standard_library
import io

standard_library.install_aliases()
from builtins import range


sp_utilities.disable_bdb_cache()




MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER = 10
# NORMALIZED_AREA_THRESHOLD_FOR_OUTLIER_DETECTION = 0.2
PERCENT_THRESHOLD_X = 0.8
PERCENT_THRESHOLD_Y = 0.2
ANGLE_ERROR_THRESHOLD = 24.0
TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND = -100
MUST_END_PROGRAM_THIS_ITERATION = -101
EMPTY_VIPER_RUN_INDICES_LIST = -102
DUMMY_INDEX_USED_AS_BUFFER = -103

NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "main"
NAME_OF_PARAMS_FILE = "rotated_reduced_params.txt"
DIR_DELIM = os.sep
AVG_VOL_FILE = "average_volume.hdf"
VAR_VOL_FILE = "variance_volume.hdf"


def calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(
    no_of_viper_runs_analyzed_together,
    no_of_viper_runs_analyzed_together_from_user_options,
    masterdir,
    rviper_iter,
    criterion_name,
    symc,
    runs_iter,
):

    # generate all possible combinations of (no_of_viper_runs_analyzed_together - 1) taken (3 - 1) at a time

    number_of_additional_combinations_for_this_viper_iteration = sp_utilities.combinations_of_n_taken_by_k(
        no_of_viper_runs_analyzed_together - 1,
        no_of_viper_runs_analyzed_together_from_user_options - 1,
    )

    criterion_measure = [
        0.0
    ] * number_of_additional_combinations_for_this_viper_iteration
    all_n_minus_1_combinations_taken_k_minus_1_at_a_time = list(
        itertools.combinations(
            list(range(no_of_viper_runs_analyzed_together - 1)),
            no_of_viper_runs_analyzed_together_from_user_options - 1,
        )
    )

    no_of_processors = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)
    my_rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)

    for idx, tuple_of_projection_indices in enumerate(
        all_n_minus_1_combinations_taken_k_minus_1_at_a_time
    ):
        if my_rank == idx % no_of_processors:
            list_of_viper_run_indices = list(tuple_of_projection_indices) + [
                no_of_viper_runs_analyzed_together - 1
            ]
            criterion_measure[idx] = measure_for_outlier_criterion(
                criterion_name, masterdir, rviper_iter, list_of_viper_run_indices, symc
            )
            plot_errors_between_any_number_of_projections(
                masterdir,
                rviper_iter,
                list_of_viper_run_indices,
                criterion_measure[idx],
                symc,
                idx,
                runs_iter,
            )

    criterion_measure = mpi.mpi_reduce(
        criterion_measure,
        number_of_additional_combinations_for_this_viper_iteration,
        mpi.MPI_FLOAT,
        mpi.MPI_SUM,
        0,
        mpi.MPI_COMM_WORLD,
    )

    if my_rank == 0:
        index_of_sorted_criterion_measure_list = [
            i[0]
            for i in sorted(
                enumerate(criterion_measure), reverse=False, key=lambda x: x[1]
            )
        ]

        list_of_viper_run_indices_for_the_current_rrr_viper_iteration = list(
            all_n_minus_1_combinations_taken_k_minus_1_at_a_time[
                index_of_sorted_criterion_measure_list[0]
            ]
        ) + [no_of_viper_runs_analyzed_together - 1]

        mainoutputdir = (
            masterdir
            + DIR_DELIM
            + NAME_OF_MAIN_DIR
            + ("%03d" + DIR_DELIM) % (rviper_iter)
        )

        if (
            criterion_measure[index_of_sorted_criterion_measure_list[0]]
            == TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND
        ):
            list_of_viper_run_indices_for_the_current_rrr_viper_iteration.insert(
                0, MUST_END_PROGRAM_THIS_ITERATION
            )
        else:
            list_of_viper_run_indices_for_the_current_rrr_viper_iteration.insert(
                0, DUMMY_INDEX_USED_AS_BUFFER
            )
            if criterion_name == "80th percentile":
                pass_criterion = (
                    criterion_measure[index_of_sorted_criterion_measure_list[0]]
                    < PERCENT_THRESHOLD_Y
                )
            elif criterion_name == "fastest increase in the last quartile":
                pass_criterion = (
                    criterion_measure[index_of_sorted_criterion_measure_list[-1]]
                    > PERCENT_THRESHOLD_Y
                )
            else:
                pass_criterion = False

            if not pass_criterion:
                list_of_viper_run_indices_for_the_current_rrr_viper_iteration = [
                    EMPTY_VIPER_RUN_INDICES_LIST
                ]

        f = open(
            mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json",
            "w",
        )
        # workaround to fix the problem of json.dump due to non-recognition of numpy.int64 objects.
        casted_list_of_viper_run_indices_for_the_current_rrr_viper_iteration = [int(v) for v in list_of_viper_run_indices_for_the_current_rrr_viper_iteration]
        json.dump(casted_list_of_viper_run_indices_for_the_current_rrr_viper_iteration[1:], f)
        f.close()

        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        return list_of_viper_run_indices_for_the_current_rrr_viper_iteration

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    return [EMPTY_VIPER_RUN_INDICES_LIST]


def identify_outliers(
    myid,
    main_node,
    rviper_iter,
    no_of_viper_runs_analyzed_together,
    no_of_viper_runs_analyzed_together_from_user_options,
    masterdir,
    bdb_stack_location,
    outlier_percentile,
    criterion_name,
    outlier_index_threshold_method,
    angle_threshold,
    symc,
    options,
    runs_iter,
):

    no_of_viper_runs_analyzed_together_must_be_incremented = 0
    do_calculation = 1

    if myid == main_node:
        mainoutputdir = (
            masterdir
            + DIR_DELIM
            + NAME_OF_MAIN_DIR
            + ("%03d" + DIR_DELIM) % (rviper_iter)
        )
        if os.path.exists(
            mainoutputdir
            + DIR_DELIM
            + "list_of_viper_runs_included_in_outlier_elimination.json"
        ):
            # list_of_independent_viper_run_indices_used_for_outlier_elimination = map(int, read_text_file(mainoutputdir + DIR_DELIM + "list_of_viper_runs_included_in_outlier_elimination.txt"))
            f = open(
                mainoutputdir
                + "list_of_viper_runs_included_in_outlier_elimination.json",
                "r",
            )
            list_of_independent_viper_run_indices_used_for_outlier_elimination = json.load(
                f
            )
            f.close()
            do_calculation = 0
        do_calculation = mpi.mpi_bcast(
            do_calculation, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD
        )[0]
    else:
        do_calculation = mpi.mpi_bcast(
            do_calculation, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD
        )[0]

    if do_calculation:
        list_of_independent_viper_run_indices_used_for_outlier_elimination = calculate_list_of_independent_viper_run_indices_used_for_outlier_elimination(
            no_of_viper_runs_analyzed_together,
            no_of_viper_runs_analyzed_together_from_user_options,
            masterdir,
            rviper_iter,
            criterion_name,
            symc,
            runs_iter,
        )

    # only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination
    # only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination
    # only master has the actual list: list_of_independent_viper_run_indices_used_for_outlier_elimination

    error_status = 0
    if myid == main_node:
        # if len(list_of_independent_viper_run_indices_used_for_outlier_elimination) == 0:
        if (
            list_of_independent_viper_run_indices_used_for_outlier_elimination[0]
            == EMPTY_VIPER_RUN_INDICES_LIST
        ):
            if (
                no_of_viper_runs_analyzed_together
                > MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER
            ):
                error_status = 1
                sp_global_def.sxprint(
                    "RVIPER reached maximum number of VIPER runs analyzed together without finding a core set of stable projections for the current RVIPER iteration (%d)! Finishing."
                    % rviper_iter
                )
                cmd = "{} {}".format(
                    "mkdir ",
                    masterdir + "MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER__Reached",
                )
                junk = sp_utilities.cmdexecute(cmd)
            else:
                # No set of solutions has been found to make a selection for outlier elimination.
                # A new independent viper run will be performed
                no_of_viper_runs_analyzed_together_must_be_incremented = 1
                cmd = "{} {}".format(
                    "rm ",
                    mainoutputdir
                    + "list_of_viper_runs_included_in_outlier_elimination.json",
                )
                junk = sp_utilities.cmdexecute(cmd)

        else:
            # Outliers are eliminated based on the viper runs contained in "list_of_independent_viper_run_indices_used_for_outlier_elimination"
            if (
                list_of_independent_viper_run_indices_used_for_outlier_elimination[0]
                == MUST_END_PROGRAM_THIS_ITERATION
            ):
                no_of_viper_runs_analyzed_together_must_be_incremented = (
                    MUST_END_PROGRAM_THIS_ITERATION
                )
                found_outliers(
                    list_of_independent_viper_run_indices_used_for_outlier_elimination[
                        1:
                    ],
                    outlier_percentile,
                    rviper_iter,
                    masterdir,
                    bdb_stack_location,
                    "use all images",
                    angle_threshold,
                    symc,
                    options,
                    runs_iter,
                )
            else:
                # still need to eliminate DUMMY_INDEX_USED_AS_BUFFER
                found_outliers(
                    list_of_independent_viper_run_indices_used_for_outlier_elimination[
                        1:
                    ],
                    outlier_percentile,
                    rviper_iter,
                    masterdir,
                    bdb_stack_location,
                    outlier_index_threshold_method,
                    angle_threshold,
                    symc,
                    options,
                    runs_iter,
                )

    sp_utilities.if_error_then_all_processes_exit_program(error_status)

    no_of_viper_runs_analyzed_together_must_be_incremented = mpi.mpi_bcast(
        no_of_viper_runs_analyzed_together_must_be_incremented,
        1,
        mpi.MPI_INT,
        0,
        mpi.MPI_COMM_WORLD,
    )[0]

    return no_of_viper_runs_analyzed_together_must_be_incremented


def plot_errors_between_any_number_of_projections(
    masterdir,
    rviper_iter,
    list_of_projection_indices,
    error_value,
    symc,
    index,
    runs_iter,
):
    matplotlib.pyplot.matplotlib.use("agg")

    # # for debugging purposes
    # if "counter" not in plot_errors_between_any_number_of_projections.__dict__:
    # 	plot_errors_between_any_number_of_projections.counter = 0
    # plot_errors_between_any_number_of_projections.counter += 1

    # main_iterations = [NAME_OF_MAIN_DIR + "%03d" % i for i in range(1, rviper_iter + 1)]
    main_iterations = [
        NAME_OF_MAIN_DIR + "%03d" % i for i in range(rviper_iter, rviper_iter + 1)
    ]

    mainoutputdir = masterdir + DIR_DELIM + main_iterations[0] + DIR_DELIM
    p = []
    for i1 in list_of_projection_indices:
        p.append(
            sp_utilities.read_text_row(
                mainoutputdir
                + NAME_OF_RUN_DIR
                + "%03d" % (i1)
                + DIR_DELIM
                + "params.txt"
            )
        )

    ti1, ti3, out = sp_multi_shc.find_common_subset(p, 0, symmetry_class=symc)
    u = []
    for i in range(len(ti3)):
        u.append([ti3[i], i])
    u.sort()
    # EMAN2.display([range(len(u)),[u[i][0] for i in xrange(len(u))]])
    matplotlib.pyplot.plot(list(range(len(u))), [u[i][0] for i in range(len(u))])

    # import json; f = open("error_curve%03d.json"%plot_errors_between_any_number_of_projections.counter, 'w')
    # json.dump([u[i][0] for i in xrange(len(u))],f); f.close()

    matplotlib.pyplot.ylabel("Error")
    matplotlib.pyplot.xlabel("Image index")
    matplotlib.pyplot.title(
        "Sorted errors between projections\nError value = %s" % error_value
    )
    which_projections = io.StringIO()
    which_projections.write(six.u("_" + "%.6f" % error_value))
    for p_i in list_of_projection_indices:
        which_projections.write(six.u("_" + "%03d" % p_i))
    for p_i in list_of_projection_indices:
        which_projections.write(
            six.u("___" + "%03d" % get_already_processed_viper_runs.r_permutation[p_i])
        )

    prj_value = which_projections.getvalue()
    ####plt.savefig(mainoutputdir + '/sorted_errors_between_projections' + prj_value + '.png')
    ####plt.savefig(mainoutputdir + DIR_DELIM + 'sorted_errors_between_projections' + "_%03d"%(runs_iter) + "_%03d"%(index) + '.png')
    matplotlib.pyplot.savefig(
        os.path.join(
            mainoutputdir,
            "sorted_errors_between_projections"
            + "_%03d" % (runs_iter)
            + "_%03d" % (index)
            + ".png",
        )
    )
    ####print('plot_errors_between_any_number_of_projections: rviper_iter %s, masterdir %s, error_value %s, getvalue %s index %s, runs_iter %s' % ( rviper_iter, masterdir, error_value, prj_value, index, runs_iter ) )
    which_projections.close()
    matplotlib.pyplot.close()


def find_index_of_discontinuity_in_derivative(
    error_curve_func,
    list_of_projection_indices,
    mainoutputdir,
    outlier_percentile,
    runs_iter,
):

    matplotlib.pyplot.matplotlib.use("agg")

    resolution = 100
    split_point_resolution = 29
    degree_of_the_fitting_polynomial = 1

    data_set_length = len(error_curve_func)

    # split_point_set = np.linspace(0.71,.99,split_point_resolution)
    # split_point_set = np.linspace(0.71,outlier_percentile/100.0,split_point_resolution)
    split_point_set = matplotlib.pyplot.np.linspace(
        0.8, old_div(outlier_percentile, 100.0), split_point_resolution
    )

    minimum_goodness_of_fit_for_both_lines = 1e20
    optimized_split_point = -1
    for split_index, split_point in enumerate(split_point_set):
        first_line_x = list(
            map(
                int,
                matplotlib.pyplot.np.linspace(0, split_point, resolution)
                * data_set_length,
            )
        )
        first_line_y = matplotlib.pyplot.np.array(
            [error_curve_func[x] for x in first_line_x]
        )
        first_line_z = matplotlib.pyplot.np.poly1d(
            matplotlib.pyplot.np.polyfit(
                first_line_x, first_line_y, degree_of_the_fitting_polynomial
            )
        )

        second_line_x = list(
            map(
                int,
                matplotlib.pyplot.np.linspace(split_point, 1, resolution)
                * data_set_length,
            )
        )
        second_line_y = matplotlib.pyplot.np.array(
            [error_curve_func[x - 1] for x in second_line_x]
        )
        second_line_z = matplotlib.pyplot.np.poly1d(
            matplotlib.pyplot.np.polyfit(
                second_line_x, second_line_y, degree_of_the_fitting_polynomial
            )
        )

        goodness_of_fit_for_both_lines = sum(
            (first_line_z(first_line_x) - first_line_y) ** 2
        )
        goodness_of_fit_for_both_lines += sum(
            (second_line_z(second_line_x) - second_line_y) ** 2
        )
        # goodness_of_fit_for_both_lines = angle((1,first_line_z[1]), (1,second_line_z[1]))
        if goodness_of_fit_for_both_lines < minimum_goodness_of_fit_for_both_lines:
            minimum_goodness_of_fit_for_both_lines = goodness_of_fit_for_both_lines
            optimized_split_point = split_point

        # split_point = optimized_split_point
        matplotlib.pyplot.plot(list(range(len(error_curve_func))), error_curve_func)

        first_line_x = list(
            map(
                int,
                matplotlib.pyplot.np.linspace(0, split_point, resolution)
                * data_set_length,
            )
        )
        first_line_y = matplotlib.pyplot.np.array(
            [error_curve_func[x] for x in first_line_x]
        )
        first_line_z = matplotlib.pyplot.np.poly1d(
            matplotlib.pyplot.np.polyfit(
                first_line_x, first_line_y, degree_of_the_fitting_polynomial
            )
        )

        matplotlib.pyplot.plot(first_line_x, first_line_z(first_line_x))

        second_line_x = list(
            map(
                int,
                matplotlib.pyplot.np.linspace(split_point, 1, resolution)
                * data_set_length,
            )
        )
        second_line_y = matplotlib.pyplot.np.array(
            [error_curve_func[x - 1] for x in second_line_x]
        )
        second_line_z = matplotlib.pyplot.np.poly1d(
            matplotlib.pyplot.np.polyfit(
                second_line_x, second_line_y, degree_of_the_fitting_polynomial
            )
        )
        matplotlib.pyplot.plot(second_line_x, second_line_z(second_line_x))

        which_projections = io.StringIO()
        which_projections.write(
            six.u("_" + "%.03f__%.6f" % (split_point, goodness_of_fit_for_both_lines))
        )
        for p_i in list_of_projection_indices:
            which_projections.write(six.u("_" + "%03d" % p_i))
        for p_i in list_of_projection_indices:
            which_projections.write(
                six.u(
                    "___" + "%03d" % get_already_processed_viper_runs.r_permutation[p_i]
                )
            )

        matplotlib.pyplot.title(
            "Sorted errors\nSplit point = %.03f\nGoodness of fit = %.6f"
            % (split_point, goodness_of_fit_for_both_lines)
        )
        ####plt.title(mainoutputdir + '/sorted_errors' + which_projections.getvalue() + '.png')
        plotfile = os.path.join(
            mainoutputdir,
            "sorted_errors" + "_%03d" % (runs_iter) + "_%03d" % (split_index) + ".png",
        )
        matplotlib.pyplot.savefig(plotfile)
        matplotlib.pyplot.close()

    split_point = optimized_split_point
    matplotlib.pyplot.plot(list(range(len(error_curve_func))), error_curve_func)

    first_line_x = list(
        map(
            int,
            matplotlib.pyplot.np.linspace(0, split_point, resolution) * data_set_length,
        )
    )
    first_line_y = matplotlib.pyplot.np.array(
        [error_curve_func[x] for x in first_line_x]
    )
    first_line_z = matplotlib.pyplot.np.poly1d(
        matplotlib.pyplot.np.polyfit(
            first_line_x, first_line_y, degree_of_the_fitting_polynomial
        )
    )

    matplotlib.pyplot.plot(first_line_x, first_line_z(first_line_x))

    second_line_x = list(
        map(
            int,
            matplotlib.pyplot.np.linspace(split_point, 1, resolution) * data_set_length,
        )
    )
    second_line_y = matplotlib.pyplot.np.array(
        [error_curve_func[x - 1] for x in second_line_x]
    )
    second_line_z = matplotlib.pyplot.np.poly1d(
        matplotlib.pyplot.np.polyfit(
            second_line_x, second_line_y, degree_of_the_fitting_polynomial
        )
    )
    matplotlib.pyplot.plot(second_line_x, second_line_z(second_line_x))

    which_projections = io.StringIO()
    which_projections.write(six.u("_" + "%.03f" % split_point))
    for p_i in list_of_projection_indices:
        which_projections.write(six.u("_" + "%03d" % p_i))
    for p_i in list_of_projection_indices:
        which_projections.write(
            six.u("___" + "%03d" % get_already_processed_viper_runs.r_permutation[p_i])
        )

    ####plt.title(mainoutputdir + '/optimized_errors' + which_projections.getvalue() + '.png')
    matplotlib.pyplot.title("Optimized error = %s" % split_point)
    plotfile = os.path.join(
        mainoutputdir, "optimized_errors" + "_%03d" % (runs_iter) + ".png"
    )
    matplotlib.pyplot.savefig(plotfile)
    matplotlib.pyplot.close()

    if optimized_split_point < 0:
        return -1

    return int(optimized_split_point * data_set_length)


def measure_for_outlier_criterion(
    criterion_name, masterdir, rviper_iter, list_of_viper_run_indices, symc
):

    # main_iterations = [NAME_OF_MAIN_DIR + "%03d" % i for i in range(1, rviper_iter + 1)]
    main_iterations = [
        NAME_OF_MAIN_DIR + "%03d" % i for i in range(rviper_iter, rviper_iter + 1)
    ]
    mainoutputdir = masterdir + DIR_DELIM + main_iterations[0] + DIR_DELIM

    p = []
    for i1 in list_of_viper_run_indices:
        p.append(
            sp_utilities.read_text_row(
                mainoutputdir
                + NAME_OF_RUN_DIR
                + "%03d" % (i1)
                + DIR_DELIM
                + "params.txt"
            )
        )
    subset, avg_diff_per_image, outp = sp_multi_shc.find_common_subset(
        p, 0, symmetry_class=symc
    )

    avg_diff_per_image.sort()
    x1 = len(avg_diff_per_image)
    y1 = avg_diff_per_image[-1]

    if y1 <= ANGLE_ERROR_THRESHOLD:
        return TRIPLET_WITH_ANGLE_ERROR_LESS_THAN_THRESHOLD_HAS_BEEN_FOUND

    if criterion_name == "80th percentile":
        return old_div(avg_diff_per_image[int(x1 * PERCENT_THRESHOLD_X)], y1)
    elif criterion_name == "fastest increase in the last quartile":
        for k in range(5, 6):
            avg_diff_per_image_diff = [
                x - avg_diff_per_image[i - k] for i, x in enumerate(avg_diff_per_image)
            ][k:]

            avg_diff_per_image_diff_max = max(avg_diff_per_image_diff)
            avg_diff_per_image_diff_max_normalized = old_div(
                max(avg_diff_per_image_diff), y1
            )

            if avg_diff_per_image_diff.index(avg_diff_per_image_diff_max) >= int(
                x1 * 0.75
            ):
                return avg_diff_per_image_diff_max_normalized
            return 0.0
    else:
        sp_global_def.sxprint("Error, no criterion name is specified!")
        mpi.mpi_finalize()
        sys.exit()


def found_outliers(
    list_of_projection_indices,
    outlier_percentile,
    rviper_iter,
    masterdir,
    bdb_stack_location,
    outlier_index_threshold_method,
    angle_threshold,
    symc,
    options,
    runs_iter,
):

    # sxheader.py bdb:nj  --consecutive  --params=OID

    mainoutputdir = (
        masterdir + DIR_DELIM + NAME_OF_MAIN_DIR + ("%03d" + DIR_DELIM) % (rviper_iter)
    )

    # if this data analysis step was already performed in the past then return
    for check_run in list_of_projection_indices:
        if not (
            os.path.exists(
                mainoutputdir
                + DIR_DELIM
                + NAME_OF_RUN_DIR
                + "%03d" % (check_run)
                + DIR_DELIM
                + NAME_OF_PARAMS_FILE
            )
        ):
            break
    else:
        return

    projs = []
    for i1 in list_of_projection_indices:
        projs.append(
            sp_utilities.read_text_row(
                mainoutputdir
                + NAME_OF_RUN_DIR
                + "%03d" % (i1)
                + DIR_DELIM
                + "params.txt"
            )
        )

    subset, avg_diff_per_image, rotated_params = sp_multi_shc.find_common_subset(
        projs, target_threshold=0, symmetry_class=symc
    )

    error_values_and_indices = []
    for i in range(len(avg_diff_per_image)):
        error_values_and_indices.append([avg_diff_per_image[i], i])
    del subset, avg_diff_per_image

    error_values_and_indices.sort()

    if outlier_index_threshold_method == "discontinuity_in_derivative":
        outlier_index_threshold = find_index_of_discontinuity_in_derivative(
            [i[0] for i in error_values_and_indices],
            list_of_projection_indices,
            mainoutputdir,
            outlier_percentile,
            runs_iter,
        )
    elif outlier_index_threshold_method == "percentile":
        outlier_index_threshold = old_div(
            outlier_percentile * (len(error_values_and_indices) - 1), 100.0
        )
    elif outlier_index_threshold_method == "angle_measure":
        error_values = [i[0] for i in error_values_and_indices]
        outlier_index_threshold = min(
            list(range(len(error_values))),
            key=lambda i: abs(error_values[i] - angle_threshold),
        )
    elif outlier_index_threshold_method == "use all images":
        outlier_index_threshold = len(error_values_and_indices)

    index_keep_images = [
        i[1] for i in error_values_and_indices[:outlier_index_threshold]
    ]
    index_outliers = [i[1] for i in error_values_and_indices[outlier_index_threshold:]]

    # print "error_values_and_indices: %f"%error_values_and_indices
    sp_global_def.sxprint(
        "outlier_index_threshold_method: ", outlier_index_threshold_method
    )
    sp_global_def.sxprint("error_values_and_indices: ", error_values_and_indices)
    sp_global_def.sxprint("index_outliers: ", index_outliers)

    reversed_sorted_index_outliers = copy.deepcopy(index_outliers)
    reversed_sorted_index_outliers.sort(reverse=True)

    for k in range(len(projs)):
        for l in reversed_sorted_index_outliers:
            del rotated_params[k][l]

    index_outliers.sort()
    index_keep_images.sort()

    sp_utilities.write_text_file(
        index_outliers, mainoutputdir + "this_iteration_index_outliers.txt"
    )
    sp_utilities.write_text_file(
        index_keep_images, mainoutputdir + "this_iteration_index_keep_images.txt"
    )

    # if len(index_outliers) < 3:
    # return False

    if len(index_outliers) > 0:
        cmd = "{} {} {} {}".format(
            "e2bdb.py ",
            bdb_stack_location + "_%03d" % (rviper_iter - 1),
            "--makevstack=" + bdb_stack_location + "_outliers_%03d" % (rviper_iter),
            "--list=" + mainoutputdir + "this_iteration_index_outliers.txt",
        )
        junk = sp_utilities.cmdexecute(cmd)
    cmd = "{} {} {} {}".format(
        "e2bdb.py ",
        bdb_stack_location + "_%03d" % (rviper_iter - 1),
        "--makevstack=" + bdb_stack_location + "_%03d" % (rviper_iter),
        "--list=" + mainoutputdir + "this_iteration_index_keep_images.txt",
    )
    junk = sp_utilities.cmdexecute(cmd)
    dat = EMAN2_cppwrap.EMData.read_images(
        bdb_stack_location + "_%03d" % (rviper_iter - 1)
    )

    sp_utilities.write_text_file(
        [dat[i].get_attr("original_image_index") for i in index_outliers],
        mainoutputdir + "index_outliers.txt",
    )
    sp_utilities.write_text_file(
        [dat[i].get_attr("original_image_index") for i in index_keep_images],
        mainoutputdir + "index_keep_images.txt",
    )

    sp_global_def.sxprint("index_outliers:: " + str(index_outliers))

    for i1 in range(len(list_of_projection_indices)):
        independent_run_dir = (
            mainoutputdir + NAME_OF_RUN_DIR + "%03d" % (list_of_projection_indices[i1])
        )
        sp_utilities.write_text_row(
            rotated_params[i1],
            mainoutputdir
            + NAME_OF_RUN_DIR
            + "%03d" % (list_of_projection_indices[i1])
            + DIR_DELIM
            + NAME_OF_PARAMS_FILE,
        )  # 5 columns

    return True


def calculate_volumes_after_rotation_and_save_them(
    ali3d_options,
    rviper_iter,
    masterdir,
    bdb_stack_location,
    mpi_rank,
    mpi_size,
    no_of_viper_runs_analyzed_together,
    no_of_viper_runs_analyzed_together_from_user_options,
    mpi_comm=-1,
):

    # This function takes into account the case in which there are more processors than images

    if mpi_comm == -1:
        mpi_comm = mpi.MPI_COMM_WORLD

    # some arguments are for debugging purposes

    mainoutputdir = (
        masterdir + DIR_DELIM + NAME_OF_MAIN_DIR + ("%03d" + DIR_DELIM) % (rviper_iter)
    )

    # list_of_projection_indices_used_for_outlier_elimination = map(int, read_text_file(mainoutputdir + DIR_DELIM + "list_of_viper_runs_included_in_outlier_elimination.txt"))
    f = open(
        mainoutputdir + "list_of_viper_runs_included_in_outlier_elimination.json", "r"
    )
    list_of_independent_viper_run_indices_used_for_outlier_elimination = json.load(f)
    f.close()

    if len(list_of_independent_viper_run_indices_used_for_outlier_elimination) == 0:
        sp_global_def.sxprint(
            "Error: len(list_of_independent_viper_run_indices_used_for_outlier_elimination)==0"
        )
        mpi.mpi_finalize()
        sys.exit()

    # if this data analysis step was already performed in the past then return
    # for future changes make sure that the file checked is the last one to be processed !!!

    # if(os.path.exists(mainoutputdir + DIR_DELIM + NAME_OF_RUN_DIR + "%03d"%(no_of_viper_runs_analyzed_together - 1) + DIR_DELIM + "rotated_volume.hdf")):
    # check_last_run = max(get_latest_directory_increment_value(mainoutputdir, NAME_OF_RUN_DIR, start_value=0), no_of_viper_runs_analyzed_together_from_user_options)
    # if(os.path.exists(mainoutputdir + DIR_DELIM + NAME_OF_RUN_DIR + "%03d"%(check_last_run) + DIR_DELIM + "rotated_volume.hdf")):
    # 	return

    # if this data analysis step was already performed in the past then return
    for check_run in list_of_independent_viper_run_indices_used_for_outlier_elimination:
        if not (
            os.path.exists(
                mainoutputdir
                + DIR_DELIM
                + NAME_OF_RUN_DIR
                + "%03d" % (check_run)
                + DIR_DELIM
                + "rotated_volume.hdf"
            )
        ):
            break
    else:
        return

    partstack = []
    # for i1 in range(0,no_of_viper_runs_analyzed_together):
    for i1 in list_of_independent_viper_run_indices_used_for_outlier_elimination:
        partstack.append(
            mainoutputdir
            + NAME_OF_RUN_DIR
            + "%03d" % (i1)
            + DIR_DELIM
            + NAME_OF_PARAMS_FILE
        )
    partids_file_name = mainoutputdir + "this_iteration_index_keep_images.txt"

    lpartids = list(map(int, sp_utilities.read_text_file(partids_file_name)))
    n_projs = len(lpartids)

    if mpi_size > n_projs:
        # if there are more processors than images
        working = int(not (mpi_rank < n_projs))
        mpi_subcomm = mpi.mpi_comm_split(
            mpi_comm, working, mpi_rank - working * n_projs
        )
        mpi_subsize = mpi.mpi_comm_size(mpi_subcomm)
        mpi_subrank = mpi.mpi_comm_rank(mpi_subcomm)
        if mpi_rank < n_projs:

            # for i in xrange(no_of_viper_runs_analyzed_together):
            for idx, i in enumerate(
                list_of_independent_viper_run_indices_used_for_outlier_elimination
            ):
                projdata = sp_utilities.getindexdata(
                    bdb_stack_location + "_%03d" % (rviper_iter - 1),
                    partids_file_name,
                    partstack[idx],
                    mpi_rank,
                    mpi_subsize,
                )
                vol = sp_multi_shc.do_volume(
                    projdata, ali3d_options, 0, mpi_comm=mpi_subcomm
                )
                del projdata
                if mpi_rank == 0:
                    vol.write_image(
                        mainoutputdir
                        + DIR_DELIM
                        + NAME_OF_RUN_DIR
                        + "%03d" % (i)
                        + DIR_DELIM
                        + "rotated_volume.hdf"
                    )
                    # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " => "
                    line = ""
                    sp_global_def.sxprint(
                        line + "Generated rec_ref_volume_run #%01d \n" % i
                    )
                del vol

        mpi.mpi_barrier(mpi_comm)
    else:
        for idx, i in enumerate(
            list_of_independent_viper_run_indices_used_for_outlier_elimination
        ):
            projdata = sp_utilities.getindexdata(
                bdb_stack_location + "_%03d" % (rviper_iter - 1),
                partids_file_name,
                partstack[idx],
                mpi_rank,
                mpi_size,
            )
            vol = sp_multi_shc.do_volume(projdata, ali3d_options, 0, mpi_comm=mpi_comm)
            del projdata
            if mpi_rank == 0:
                vol.write_image(
                    mainoutputdir
                    + DIR_DELIM
                    + NAME_OF_RUN_DIR
                    + "%03d" % (i)
                    + DIR_DELIM
                    + "rotated_volume.hdf"
                )
                # line = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " => "
                line = ""
                sp_global_def.sxprint(line + "Generated rec_ref_volume_run #%01d" % i)
            del vol

    if mpi_rank == 0:
        # Align all rotated volumes, calculate their average and save as an overall result
        # vls = [None]*no_of_viper_runs_analyzed_together
        vls = [None] * len(
            list_of_independent_viper_run_indices_used_for_outlier_elimination
        )
        # for i in xrange(no_of_viper_runs_analyzed_together):
        for idx, i in enumerate(
            list_of_independent_viper_run_indices_used_for_outlier_elimination
        ):
            vls[idx] = sp_utilities.get_im(
                mainoutputdir
                + DIR_DELIM
                + NAME_OF_RUN_DIR
                + "%03d" % (i)
                + DIR_DELIM
                + "rotated_volume.hdf"
            )
            sp_utilities.set_params3D(vls[idx], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0])
        asa, sas = sp_statistics.ave_var(vls)
        # do the alignment
        nx = asa.get_xsize()
        radius = old_div(nx, 2) - 0.5
        st = EMAN2_cppwrap.Util.infomask(
            asa * asa, sp_utilities.model_circle(radius, nx, nx, nx), True
        )
        goal = st[0]
        going = True
        while going:
            sp_utilities.set_params3D(asa, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 1.0])
            # for i in xrange(no_of_viper_runs_analyzed_together):
            for idx, i in enumerate(
                list_of_independent_viper_run_indices_used_for_outlier_elimination
            ):
                o = sp_applications.ali_vol(
                    vls[idx], asa, 7.0, 5.0, radius
                )  # range of angles and shifts, maybe should be adjusted
                p = sp_utilities.get_params3D(o)
                del o
                sp_utilities.set_params3D(vls[idx], p)
            asa, sas = sp_statistics.ave_var(vls)
            st = EMAN2_cppwrap.Util.infomask(
                asa * asa, sp_utilities.model_circle(radius, nx, nx, nx), True
            )
            if st[0] > goal:
                goal = st[0]
            else:
                going = False
        # over and out
        asa.write_image(os.path.join(mainoutputdir, AVG_VOL_FILE))
        sas.write_image(os.path.join(mainoutputdir, VAR_VOL_FILE))
    return


def get_already_processed_viper_runs(run_get_already_processed_viper_runs):

    if run_get_already_processed_viper_runs:
        location_location = "/Users/hvoicu/Analysis/rrviper/particle__PIC_ISAC_g1_clean/0001__sim_r_viper_pool_001/"
        location_location = "/Users/hvoicu/Analysis/rrviper/particle__sp_MED_isac_clean_v1/0001__sim_r_viper_pool_001/"

        if "counter" not in get_already_processed_viper_runs.__dict__:
            # function needs to be called once before being used !
            get_already_processed_viper_runs.counter = -2

            path, dirs, files = next(os.walk(location_location))
            # dirs = filter(lambda x:'run' in x, dirs)
            dirs = [x for x in dirs if matplotlib.re.search("run\d\d\d$", x)]
            get_already_processed_viper_runs.r_permutation = list(range(len(dirs)))
            random.shuffle(get_already_processed_viper_runs.r_permutation)
            sp_global_def.sxprint(str(get_already_processed_viper_runs.r_permutation))
        get_already_processed_viper_runs.counter += 1
        sp_global_def.sxprint(
            "get_already_processed_viper_runs.counter: "
            + str(get_already_processed_viper_runs.counter)
        )
        # if get_already_processed_viper_runs.counter > 9:
        if get_already_processed_viper_runs.counter > (
            MAXIMUM_NO_OF_VIPER_RUNS_ANALYZED_TOGETHER - 1
        ):
            sp_global_def.sxprint("get_already_processed_viper_runs.counter > 9")
            mpi.mpi_finalize()
            sys.exit()

        return (
            location_location
            + NAME_OF_RUN_DIR
            + "%03d"
            % get_already_processed_viper_runs.r_permutation[
                get_already_processed_viper_runs.counter
            ]
        )
    else:
        get_already_processed_viper_runs.r_permutation = [0] * 20


def run():

    main_node = 0
    mpi_comm = mpi.MPI_COMM_WORLD
    myid = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    mpi_size = mpi.mpi_comm_size(
        mpi.MPI_COMM_WORLD
    )  # Total number of processes, passed by --np option.

    progname = os.path.basename(sys.argv[0])
    usage = (
        progname
        + " stack  [output_directory] --ir=inner_radius --radius=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1  --fl --aa --ref_a=S --sym=c1"
    )
    usage += """

stack			2D images in a stack file: (default required string)
output_directory: directory name into which the output files will be written.  If it does not exist, the directory will be created.  If it does exist, the program will continue executing from where it stopped (if it did not already reach the end). The "--use_latest_master_directory" option can be used to choose the most recent directory that starts with "master".
"""

    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--radius",
        type="int",
        default=29,
        help="radius of the particle: has to be less than < int(nx/2)-1 (default 29)",
    )

    parser.add_option(
        "--ir",
        type="int",
        default=1,
        help="inner radius for rotational search: > 0 (default 1)",
    )
    parser.add_option(
        "--rs",
        type="int",
        default=1,
        help="step between rings in rotational search: >0 (default 1)",
    )
    parser.add_option(
        "--xr",
        type="string",
        default="0",
        help="range for translation search in x direction: search is +/xr in pixels (default '0')",
    )
    parser.add_option(
        "--yr",
        type="string",
        default="0",
        help="range for translation search in y direction: if omitted will be set to xr, search is +/yr in pixels (default '0')",
    )
    parser.add_option(
        "--ts",
        type="string",
        default="1.0",
        help="step size of the translation search in x-y directions: search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional (default '1.0')",
    )
    parser.add_option(
        "--delta",
        type="string",
        default="2.0",
        help="angular step of reference projections: (default '2.0')",
    )
    # parser.add_option("--an",       type="string", default= "-1",              help="angular neighborhood for local searches (phi and theta)")
    parser.add_option(
        "--center",
        type="float",
        default=-1.0,
        help="centering of 3D template: average shift method; 0: no centering; 1: center of gravity (default -1.0)",
    )
    parser.add_option(
        "--maxit1",
        type="int",
        default=400,
        help="maximum number of iterations performed for the GA part: (default 400)",
    )
    parser.add_option(
        "--maxit2",
        type="int",
        default=50,
        help="maximum number of iterations performed for the finishing up part: (default 50)",
    )
    parser.add_option(
        "--L2threshold",
        type="float",
        default=0.03,
        help="stopping criterion of GA: given as a maximum relative dispersion of volumes' L2 norms: (default 0.03)",
    )
    parser.add_option(
        "--doga",
        type="float",
        default=0.1,
        help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga: (default 0.1)",
    )
    parser.add_option(
        "--n_shc_runs",
        type="int",
        default=4,
        help="number of quasi-independent shc runs (same as '--nruns' parameter from sxviper.py): (default 4)",
    )
    parser.add_option(
        "--n_rv_runs",
        type="int",
        default=10,
        help="number of rviper iterations: (default 10)",
    )
    parser.add_option(
        "--n_v_runs",
        type="int",
        default=3,
        help="number of viper runs for each r_viper cycle: (default 3)",
    )
    parser.add_option(
        "--outlier_percentile",
        type="float",
        default=95.0,
        help="percentile above which outliers are removed every rviper iteration: (default 95.0)",
    )
    parser.add_option(
        "--iteration_start",
        type="int",
        default=0,
        help="starting iteration for rviper: 0 means go to the most recent one (default 0)",
    )
    # parser.add_option("--CTF",      action="store_true", default=False,        help="NOT IMPLEMENTED Consider CTF correction during the alignment ")
    # parser.add_option("--snr",      type="float",  default= 1.0,               help="Signal-to-Noise Ratio of the data (default 1.0)")
    parser.add_option(
        "--ref_a",
        type="string",
        default="S",
        help="method for generating the quasi-uniformly distributed projection directions: (default S)",
    )
    parser.add_option(
        "--sym",
        type="string",
        default="c1",
        help="point-group symmetry of the structure: (default c1)",
    )
    # parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--function", type="string", default="ref_ali3d", help=optparse.SUPPRESS_HELP
    )
    parser.add_option(
        "--npad",
        type="int",
        default=2,
        help="padding size for 3D reconstruction: (default 2)",
    )
    # parser.add_option("--npad", type="int",  default= 2,            help="padding size for 3D reconstruction (default 2)")

    # options introduced for the do_volume function
    parser.add_option(
        "--fl",
        type="float",
        default=0.25,
        help="cut-off frequency applied to the template volume: using a hyperbolic tangent low-pass filter (default 0.25)",
    )
    parser.add_option(
        "--aa",
        type="float",
        default=0.1,
        help="fall-off of hyperbolic tangent low-pass filter: (default 0.1)",
    )
    parser.add_option(
        "--pwreference",
        type="string",
        default="",
        help="text file with a reference power spectrum: (default none)",
    )
    parser.add_option(
        "--mask3D", type="string", default=None, help="3D mask file (default sphere)"
    )
    parser.add_option(
        "--moon_elimination",
        type="string",
        default="",
        help="elimination of disconnected pieces: two arguments: mass in KDa and pixel size in px/A separated by comma, no space (default none)",
    )

    # options introduced for the angular distribution
    parser.add_option(
        "--dpi",
        type="int",
        default=72,
        help="resolution of angular distribution BILD file, dpi (default 72)",
    )

    # options introduced for the angular distribution
    parser.add_option(
        "--theta1",
        type="float",
        default=-1.0,
        help="starting restriction angle for out-of-plane tilt",
    )
    parser.add_option(
        "--theta2",
        type="float",
        default=-1.0,
        help="ending restriction angle for out-of-plane tilt",
    )

    # used for debugging, help is supressed with SUPPRESS_HELP
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--my_random_seed", type="int", default=123, help=optparse.SUPPRESS_HELP
    )
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--run_get_already_processed_viper_runs",
        action="store_true",
        dest="run_get_already_processed_viper_runs",
        default=False,
        help=optparse.SUPPRESS_HELP,
    )
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--use_latest_master_directory",
        action="store_true",
        dest="use_latest_master_directory",
        default=False,
        help=optparse.SUPPRESS_HELP,
    )

    parser.add_option(
        "--criterion_name",
        type="string",
        default="80th percentile",
        help="criterion deciding if volumes have a core set of stable projections: '80th percentile', other options:'fastest increase in the last quartile' (default '80th percentile')",
    )
    parser.add_option(
        "--outlier_index_threshold_method",
        type="string",
        default="discontinuity_in_derivative",
        help="method that decides which images to keep: discontinuity_in_derivative, other options:percentile, angle_measure (default discontinuity_in_derivative)",
    )
    parser.add_option(
        "--angle_threshold",
        type="int",
        default=30,
        help="angle threshold for projection removal if using 'angle_measure': (default 30)",
    )
    parser.add_option(
        "--filament_width",
        type="int",
        default=-1,
        help="When this is set to a non-default value helical data is assumed in which case the input images will be aligned to a mask with this width. (Default: -1)",
    )
    parser.add_option(
        "--resample_ratio",
        type="string",
        default="1.0",
        help="Specify a value larger than 0.0. By default, the program does not resample the input map (i.e. resample ratio is 1.0). Use this option maily to restore the original dimensions or pixel size of VIPER or R-VIPER model. Alternatively, specify the path to the output directory of an ISAC2 run. The program automatically extracts the resampling ratio used by the ISAC2 run. (default '1.0')",
    )

    required_option_list = ["radius"]
    (options, args) = parser.parse_args(sys.argv[1:])

    try:
        float(options.resample_ratio)
    except ValueError:
        isac_shrink_file = open(
            os.path.join(options.resample_ratio, "README_shrink_ratio.txt"), "r"
        )
        isac_shrink_lines = isac_shrink_file.readlines()
        isac_shrink_ratio = float(
            isac_shrink_lines[5]
        )  # 6th line: shrink ratio (= [target particle radius]/[particle radius]) used in the ISAC run
        isac_radius = float(
            isac_shrink_lines[6]
        )  # 7th line: particle radius at original pixel size used in the ISAC run
        isac_shrink_file.close()
        if myid == main_node:
            sp_global_def.sxprint(" ")
            sp_global_def.sxprint(
                "ISAC2 run directory path is specified with --resample_ratio option..."
            )
            sp_global_def.sxprint("Extracted parameter values")
            sp_global_def.sxprint(
                "  ISAC shrink ratio    : {}".format(isac_shrink_ratio)
            )
            sp_global_def.sxprint("  ISAC particle radius : {}".format(isac_radius))
        options.resample_ratio = isac_shrink_ratio
    else:
        options.resample_ratio = float(options.resample_ratio)
        if myid == main_node:
            if options.resample_ratio != 1.0:
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "Resample ratio {} is specified with --resample_ratio option...".format(
                        options.resample_ratio
                    )
                )
            else:
                sp_global_def.sxprint(" ")
                sp_global_def.sxprint(
                    "Resample ratio is {}. The program does not resample the input volume...".format(
                        options.resample_ratio
                    )
                )

    options.CTF = False
    options.snr = 1.0
    options.an = -1

    if options.moon_elimination == "":
        options.moon_elimination = []
    else:
        options.moon_elimination = list(map(float, options.moon_elimination.split(",")))

    # Making sure all required options appeared.
    for required_option in required_option_list:
        if not options.__dict__[required_option]:
            sp_global_def.sxprint(
                "\n ==%s== mandatory option is missing.\n" % required_option
            )
            sp_global_def.sxprint(
                "Please run '" + progname + " -h' for detailed options"
            )
            sp_global_def.ERROR(
                "Invalid number of parameters used. Please see usage information above."
            )
            return

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    if myid == main_node:
        sp_global_def.sxprint(
            "****************************************************************"
        )
        EMAN2_cppwrap.Util.version()
        sp_global_def.sxprint(
            "****************************************************************"
        )
        sys.stdout.flush()
    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    # this is just for benefiting from a user friendly parameter name
    my_random_seed = options.my_random_seed
    criterion_name = options.criterion_name
    outlier_index_threshold_method = options.outlier_index_threshold_method
    use_latest_master_directory = options.use_latest_master_directory
    iteration_start_default = options.iteration_start
    number_of_rrr_viper_runs = options.n_rv_runs
    no_of_viper_runs_analyzed_together_from_user_options = options.n_v_runs
    no_of_shc_runs_analyzed_together = options.n_shc_runs
    outlier_percentile = options.outlier_percentile
    angle_threshold = options.angle_threshold
    options.ou = options.radius  # used by multi_shc
    options.method = options.ref_a  # used by angular_distribution

    run_get_already_processed_viper_runs = options.run_get_already_processed_viper_runs
    get_already_processed_viper_runs(run_get_already_processed_viper_runs)

    random.seed(my_random_seed)

    if len(args) < 1 or len(args) > 3:
        sp_global_def.sxprint("Usage: " + usage)
        sp_global_def.sxprint("Please run '" + progname + " -h' for detailed options")
        sp_global_def.ERROR(
            "Invalid number of parameters used. Please see usage information above."
        )
        return

    # if len(args) > 2:
    # 	ref_vol = get_im(args[2])
    # else:
    ref_vol = None

    # error_status = None
    # if myid == 0:
    # 	number_of_images = EMUtil.get_image_count(args[0])
    # 	if mpi_size > number_of_images:
    # 		error_status = ('Number of processes supplied by --np in mpirun needs to be less than or equal to %d (total number of images) ' % number_of_images, getframeinfo(currentframe()))
    # if_error_then_all_processes_exit_program(error_status)

    bdb_stack_location = ""

    # If necessary, append trailing '/' to directory
    masterdir = ""
    if len(args) == 2:
        masterdir = args[1]
        if masterdir[-1] != DIR_DELIM:
            masterdir += DIR_DELIM
    elif len(args) == 1:
        if use_latest_master_directory:
            all_dirs = [d for d in os.listdir("") if os.path.isdir(d)]
            r = re.compile("^master.*$")
            all_dirs = list(filter(r.match, all_dirs))
            if len(all_dirs) > 0:
                # all_dirs = max(all_dirs, key=os.path.getctime)
                masterdir = max(all_dirs, key=os.path.getmtime)
                masterdir += DIR_DELIM

    log = sp_logger.Logger(sp_logger.BaseLogger_Files())

    error_status = 0
    if mpi_size % no_of_shc_runs_analyzed_together != 0:
        sp_global_def.ERROR(
            "Number of processes needs to be a multiple of the number of quasi-independent runs (shc) within each viper run. "
            "Total quasi-independent runs by default are 3, you can change it by specifying "
            "--n_shc_runs option (in sxviper this option is called --nruns). Also, to improve communication time it is recommended that "
            "the number of processes divided by the number of quasi-independent runs is a power "
            "of 2 (e.g. 2, 4, 8 or 16 depending on how many physical cores each node has)."
        )
        error_status = 1
    sp_utilities.if_error_then_all_processes_exit_program(error_status)

    # Create folder for all results or check if there is one created already
    if myid == main_node:
        # cmd = "{}".format("Rmycounter ccc")
        # junk = cmdexecute(cmd)

        if masterdir == "":
            timestring = time.strftime(
                "%Y_%m_%d__%H_%M_%S" + DIR_DELIM, time.localtime()
            )
            masterdir = "master" + timestring

        if not os.path.exists(masterdir):
            cmd = "{} {}".format("mkdir -p", masterdir)
            junk = sp_utilities.cmdexecute(cmd)
        sp_global_def.write_command(masterdir)

        if options.filament_width != -1:
            suffix = "_000_no_align"
        else:
            suffix = "_000"
        if ":" in args[0]:
            bdb_stack_location = (
                args[0].split(":")[0] + ":" + masterdir + args[0].split(":")[1]
            )
            org_stack_location = args[0]

            if not os.path.exists(os.path.join(masterdir, "EMAN2DB" + DIR_DELIM)):
                # cmd = "{} {}".format("cp -rp EMAN2DB", masterdir, "EMAN2DB" DIR_DELIM)
                # junk = cmdexecute(cmd)
                cmd = "{} {} {}".format(
                    "e2bdb.py",
                    org_stack_location,
                    "--makevstack=" + bdb_stack_location + suffix,
                )
                junk = sp_utilities.cmdexecute(cmd)

                try:
                    sp_applications.header(
                        bdb_stack_location + suffix,
                        params="original_image_index",
                        fprint=True,
                    )
                    sp_global_def.sxprint("Images were already indexed!")
                except KeyError:
                    sp_global_def.sxprint("Indexing images")
                    sp_applications.header(
                        bdb_stack_location + suffix,
                        params="original_image_index",
                        consecutive=True,
                    )
        else:

            filename = os.path.basename(args[0])
            bdb_stack_location = "bdb:" + masterdir + os.path.splitext(filename)[0]
            bdb_path = (
                masterdir
                + "EMAN2DB"
                + DIR_DELIM
                + os.path.splitext(filename)[0]
                + suffix
                + ".bdb"
            )

            junk = True
            if not os.path.exists(bdb_path):
                cmd = "{} {} {}".format(
                    "sp_cpy.py  ", args[0], bdb_stack_location + suffix
                )
                junk = sp_utilities.cmdexecute(cmd)
                # cpy(args[0], bdb_stack_location + "_000")  # without subprocess

                if junk:
                    try:
                        sp_applications.header(
                            bdb_stack_location + suffix,
                            params="original_image_index",
                            fprint=True,
                        )
                        sp_global_def.sxprint("Images were already indexed!")
                    except KeyError:
                        sp_global_def.sxprint("Indexing images")
                        sp_applications.header(
                            bdb_stack_location + suffix,
                            params="original_image_index",
                            consecutive=True,
                        )

        if options.filament_width != -1:
            sp_global_def.sxprint("XXXXXXXXXXXXXXXXX")
            sp_global_def.sxprint("Do helical prealignment")
            sp_global_def.sxprint("XXXXXXXXXXXXXXXXX")
            helical_stack = bdb_stack_location + "_000_helical"
            transform_stack = bdb_stack_location + "_000"
            cmd = "{} {} {}".format(
                "e2bdb.py", bdb_stack_location + suffix, "--makevstack=" + helical_stack
            )
            junk = sp_utilities.cmdexecute(cmd)
            sp_applications.header(helical_stack, params="xform.align2d", zero=True)
            first_proj = sp_utilities.get_im(helical_stack)
            mask_dim = first_proj.get_xsize()
            mask = sp_utilities.model_rotated_rectangle2D(
                radius_long=mask_dim,  # long  edge of the rectangular mask
                radius_short=old_div(
                    int(options.filament_width * options.resample_ratio + 0.5), 2
                ),  # short edge of the rectangular mask
                nx=mask_dim,
                ny=mask_dim,
                angle=90,
            )
            mask_file = "{0}/filament_mask.hdf".format(masterdir)
            mask.write_image(mask_file)
            sp_applications.mref_ali2d(
                stack=helical_stack,
                refim=mask_file,
                outdir=os.path.join(masterdir, "filament_mref"),
                maskfile=mask_file,
                ir=1,
                ou=options.radius,
                center=0,
                maxit=10,
            )
            cmd = "{} {} {}".format("sp_transform2d.py", helical_stack, transform_stack)
            junk = sp_utilities.cmdexecute(cmd)
    else:
        cmd = ""
        junk = None

    junk = sp_utilities.wrap_mpi_bcast(junk, main_node, mpi.MPI_COMM_WORLD)
    if not junk:
        mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
        sp_global_def.ERROR(
            "Command failed! {0}. Exit here!".format(cmd), myid=myid, action=0
        )
        mpi.mpi_finalize()
        exit(1)

    # send masterdir to all processes
    dir_len = len(masterdir) * int(myid == main_node)
    dir_len = mpi.mpi_bcast(dir_len, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)[0]
    masterdir = sp_utilities.wrap_mpi_bcast(
        masterdir, main_node, mpi.MPI_COMM_WORLD
    )
    masterdir = "".join(masterdir)
    if masterdir[-1] != DIR_DELIM:
        masterdir += DIR_DELIM

    sp_global_def.LOGFILE = os.path.join(masterdir, sp_global_def.LOGFILE)

    # send bdb_stack_location to all processes
    dir_len = len(bdb_stack_location) * int(myid == main_node)
    dir_len = mpi.mpi_bcast(dir_len, 1, mpi.MPI_INT, 0, mpi.MPI_COMM_WORLD)[0]
    bdb_stack_location = sp_utilities.wrap_mpi_bcast(bdb_stack_location, main_node, mpi.MPI_COMM_WORLD)
    bdb_stack_location = "".join(bdb_stack_location)

    iteration_start = sp_utilities.get_latest_directory_increment_value(
        masterdir, "main"
    )

    if iteration_start_default != 0:
        iteration_start = iteration_start_default
    if myid == main_node:
        if iteration_start < iteration_start_default:
            sp_global_def.ERROR(
                "Starting iteration provided is greater than last iteration performed. Quitting program"
            )
            error_status = 1
        if number_of_rrr_viper_runs < iteration_start:
            sp_global_def.ERROR(
                "Please provide number of rviper runs (--n_rv_runs) greater than number of iterations already performed."
            )
            error_status = 1

    sp_utilities.if_error_then_all_processes_exit_program(error_status)

    symc = sp_fundamentals.symclass(options.sym)

    # Loop through iterations
    for rviper_iter in range(iteration_start, number_of_rrr_viper_runs + 1):
        if myid == main_node:
            all_projs = EMAN2_cppwrap.EMData.read_images(
                bdb_stack_location + "_%03d" % (rviper_iter - 1)
            )
            sp_global_def.sxprint("XXXXXXXXXXXXXXXXX")
            sp_global_def.sxprint(
                "Number of projections (in loop): " + str(len(all_projs))
            )
            sp_global_def.sxprint("XXXXXXXXXXXXXXXXX")
            subset = list(range(len(all_projs)))
        else:
            all_projs = None
            subset = None

        runs_iter = (
            sp_utilities.get_latest_directory_increment_value(
                masterdir + NAME_OF_MAIN_DIR + "%03d" % rviper_iter,
                DIR_DELIM + NAME_OF_RUN_DIR,
                start_value=0,
            )
            - 1
        )
        no_of_viper_runs_analyzed_together = max(
            runs_iter + 2, no_of_viper_runs_analyzed_together_from_user_options
        )

        first_time_entering_the_loop_need_to_do_full_check_up = True
        while True:
            runs_iter += 1

            if not first_time_entering_the_loop_need_to_do_full_check_up:
                if runs_iter >= no_of_viper_runs_analyzed_together:
                    break
            first_time_entering_the_loop_need_to_do_full_check_up = False

            this_run_is_NOT_complete = 0
            if myid == main_node:
                independent_run_dir = (
                    masterdir
                    + DIR_DELIM
                    + NAME_OF_MAIN_DIR
                    + ("%03d" + DIR_DELIM + NAME_OF_RUN_DIR + "%03d" + DIR_DELIM)
                    % (rviper_iter, runs_iter)
                )
                if run_get_already_processed_viper_runs:
                    cmd = "{} {}".format(
                        "mkdir -p",
                        masterdir
                        + DIR_DELIM
                        + NAME_OF_MAIN_DIR
                        + ("%03d" + DIR_DELIM) % (rviper_iter),
                    )
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format("rm -rf", independent_run_dir)
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format(
                        "cp -r",
                        get_already_processed_viper_runs() + " " + independent_run_dir,
                    )
                    junk = sp_utilities.cmdexecute(cmd)

                if os.path.exists(independent_run_dir + "log.txt") and (
                    sp_utilities.string_found_in_file(
                        "Finish VIPER2", independent_run_dir + "log.txt"
                    )
                ):
                    this_run_is_NOT_complete = 0
                else:
                    this_run_is_NOT_complete = 1
                    cmd = "{} {}".format("rm -rf", independent_run_dir)
                    junk = sp_utilities.cmdexecute(cmd)
                    cmd = "{} {}".format("mkdir -p", independent_run_dir)
                    junk = sp_utilities.cmdexecute(cmd)

                this_run_is_NOT_complete = mpi.mpi_bcast(
                    this_run_is_NOT_complete,
                    1,
                    mpi.MPI_INT,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )[0]
                dir_len = len(independent_run_dir)
                dir_len = mpi.mpi_bcast(
                    dir_len, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD
                )[0]
                independent_run_dir = sp_utilities.wrap_mpi_bcast(
                    independent_run_dir,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )
                independent_run_dir = "".join(independent_run_dir)
            else:
                this_run_is_NOT_complete = mpi.mpi_bcast(
                    this_run_is_NOT_complete,
                    1,
                    mpi.MPI_INT,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )[0]
                dir_len = 0
                independent_run_dir = ""
                dir_len = mpi.mpi_bcast(
                    dir_len, 1, mpi.MPI_INT, main_node, mpi.MPI_COMM_WORLD
                )[0]
                independent_run_dir = sp_utilities.wrap_mpi_bcast(
                    independent_run_dir,
                    main_node,
                    mpi.MPI_COMM_WORLD,
                )
                independent_run_dir = "".join(independent_run_dir)

            if this_run_is_NOT_complete:
                mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

                if independent_run_dir[-1] != DIR_DELIM:
                    independent_run_dir += DIR_DELIM

                log.prefix = independent_run_dir

                options.user_func = sp_user_functions.factory[options.function]

                # for debugging purposes
                # if (myid == main_node):
                # cmd = "{} {}".format("cp ~/log.txt ", independent_run_dir)
                # junk = cmdexecute(cmd)
                # cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "param%03d.txt"%runs_iter)
                # cmd = "{} {}{}".format("cp ~/paramdir/params$(mycounter ccc).txt ", independent_run_dir, "params.txt")
                # junk = cmdexecute(cmd)

                if myid == main_node:
                    sp_utilities.store_value_of_simple_vars_in_json_file(
                        masterdir + "program_state_stack.json",
                        locals(),
                        exclude_list_of_vars=["usage"],
                        vars_that_will_show_only_size=["subset"],
                    )
                    sp_utilities.store_value_of_simple_vars_in_json_file(
                        masterdir + "program_state_stack.json",
                        options.__dict__,
                        write_or_append="a",
                    )

                out_params, out_vol, out_peaks = sp_multi_shc.multi_shc(
                    all_projs,
                    subset,
                    no_of_shc_runs_analyzed_together,
                    options,
                    mpi_comm=mpi_comm,
                    log=log,
                    ref_vol=ref_vol,
                )
                # options used (not all below are accepted as input from command line):
                # by multi_shc: sym, delta
                # 	-> do_volume: sym, npad, CTF, snr, mask3D or ou, pwreference, fl, aa
                # 	-> ali3d_multishc: ir, rs, ou, xr, yr, ts, an, delta, doga, center, CTF, ref_a, L2threshold, maxit1, moon_elimination
                # 	-> model_circle: ou
                # 	-> ali3d_multishc_2: ir, rs, ou, xr, yr, ts, an, delta, center, CTF, ref_a, maxit2

                # end of: if this_run_is_NOT_complete:

            if runs_iter >= (no_of_viper_runs_analyzed_together_from_user_options - 1):
                increment_for_current_iteration = identify_outliers(
                    myid,
                    main_node,
                    rviper_iter,
                    no_of_viper_runs_analyzed_together,
                    no_of_viper_runs_analyzed_together_from_user_options,
                    masterdir,
                    bdb_stack_location,
                    outlier_percentile,
                    criterion_name,
                    outlier_index_threshold_method,
                    angle_threshold,
                    symc,
                    options,
                    runs_iter,
                )

                if increment_for_current_iteration == MUST_END_PROGRAM_THIS_ITERATION:
                    break

                no_of_viper_runs_analyzed_together += increment_for_current_iteration

        # end of independent viper loop

        calculate_volumes_after_rotation_and_save_them(
            options,
            rviper_iter,
            masterdir,
            bdb_stack_location,
            myid,
            mpi_size,
            no_of_viper_runs_analyzed_together,
            no_of_viper_runs_analyzed_together_from_user_options,
        )

        if increment_for_current_iteration == MUST_END_PROGRAM_THIS_ITERATION:
            if myid == main_node:
                sp_global_def.sxprint(
                    "RVIPER found a core set of stable projections for the current RVIPER iteration (%d), the maximum angle difference between corresponding projections from different VIPER volumes is less than %.2f. Finishing."
                    % (rviper_iter, ANGLE_ERROR_THRESHOLD)
                )
                iter_avg = (
                    masterdir
                    + DIR_DELIM
                    + NAME_OF_MAIN_DIR
                    + "%03d" % (rviper_iter)
                    + DIR_DELIM
                    + AVG_VOL_FILE
                )
                iter_var = (
                    masterdir
                    + DIR_DELIM
                    + NAME_OF_MAIN_DIR
                    + "%03d" % (rviper_iter)
                    + DIR_DELIM
                    + VAR_VOL_FILE
                )
                master_avg = os.path.join(
                    masterdir, "average_volume_%03d.hdf" % rviper_iter
                )
                master_var = os.path.join(
                    masterdir, "variance_volume_%03d.hdf" % rviper_iter
                )
                sp_global_def.sxprint(
                    "Copying average and variance from iteration %03d to output directory %s"
                    % (rviper_iter, masterdir)
                )
                shutil.copyfile(iter_avg, master_avg)
                shutil.copyfile(iter_var, master_var)
            break
    else:
        if myid == main_node:
            sp_global_def.sxprint(
                "After running the last iteration (%d), RVIPER did not find a set of projections with the maximum angle difference between corresponding projections from different VIPER volumes less than %.2f Finishing."
                % (rviper_iter, ANGLE_ERROR_THRESHOLD)
            )

    # end of RVIPER iteration loop

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

def main():
    mpi.mpi_init(0, [])
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()


if __name__ == "__main__":
    main()
