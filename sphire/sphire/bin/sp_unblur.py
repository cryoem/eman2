#! /usr/bin/env python
from __future__ import print_function
from __future__ import division
from past.utils import old_div

#
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


import argparse
import glob
import mpi
import os
from ..libpy import sp_global_def
from ..libpy import sp_utilities
from ..libpy import sp_applications
import subprocess
import time

mpi.mpi_init(0, [])


def parse_args():
    """
	Parse command line arguments.

	Arguments:
	None

	Returns:
	dictionary of command line arguments
	"""

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Mandatory arguments
    parser.add_argument(
        "unblur_path", help="Path to the unblur executable that is coming with cisTEM."
    )
    parser.add_argument(
        "input_micrograph_pattern",
        help="Pattern of the micrographs to process, the variable part of the name should be replaced by the wildcard symbol * e.g. Movies/prefix_*_suffix.mrcs",
    )
    parser.add_argument("output_directory", help="Output directory for the results")

    # Optional arguments
    parser.add_argument(
        "--selection_file", default=None, help="Selection list for the micrographs."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Do not check if the output directory already exists. Enabling this option can lead to data loss.",
    )

    group_general = parser.add_argument_group("General")
    group_general.add_argument(
        "--pixel_size", default=1.0, help="Pixel size of images (A)"
    )
    group_general.add_argument("--bin_factor", default=1.0, help="Output bin factor")

    group_dose = parser.add_argument_group("Dose adjustment")
    group_dose.add_argument(
        "--skip_dose_adjustment",
        action="store_true",
        default=False,
        help="Do not apply exposure filter?",
    )
    group_dose.add_argument(
        "--additional_dose_unadjusted",
        action="store_true",
        default=False,
        help="Run unblur twice to produce undoseweighted images.",
    )
    group_dose.add_argument(
        "--voltage", default=300.0, help="Acceleration voltage (kV)"
    )
    group_dose.add_argument(
        "--exposure_per_frame", default=1.0, help="Exposure per frame (e/A^2)"
    )
    group_dose.add_argument(
        "--pre_exposure", default=0.0, help="Pre-exposure amount (e/A^2)"
    )

    group_expert = parser.add_argument_group("Expert options")
    group_expert.add_argument(
        "--min_shift_initial", default=2.0, help="Minimum shift for initial search (A)"
    )
    group_expert.add_argument(
        "--outer_radius", default=20.0, help="Outer radius shift limit (A)"
    )
    group_expert.add_argument(
        "--b_factor", default=1500.0, help="B-factor to apply to images (A^2)"
    )
    group_expert.add_argument(
        "--half_width_vert", default=1, help="Half-width of vertical Fourier mask"
    )
    group_expert.add_argument(
        "--half_width_hor", default=1, help="Half-width of horizontal Fourier mask"
    )
    group_expert.add_argument(
        "--termination", default=1, help="Termination shift threshold (A)"
    )
    group_expert.add_argument(
        "--max_iterations", default=20, help="Maximum number of iterations"
    )
    group_expert.add_argument(
        "--dont_restore_noise_power",
        action="store_true",
        default=False,
        help="Do not Restore Noise Power?",
    )
    group_expert.add_argument("--gain_file", default=None, help="Gain image filename")
    group_expert.add_argument(
        "--first_frame", default=1, help="First frame to use for sum"
    )
    group_expert.add_argument(
        "--last_frame", default=0, help="Last frame to use for sum (0 for last frame)"
    )

    group_mag = parser.add_argument_group("Magnification correction")
    group_mag.add_argument(
        "--distortion_angle", default=0.0, help="Distortion Angle (Degrees)"
    )
    group_mag.add_argument("--major_scale", default=1.0, help="Major Scale")
    group_mag.add_argument("--minor_scale", default=1.0, help="Minor Scale")

    args = vars(parser.parse_args())

    check_groups = ("Expert options", "Magnification correction")
    for entry in check_groups:
        args[entry] = False

    for group in parser._action_groups:
        if group.title in check_groups:
            for argument in group._group_actions:
                if argument.default != args[argument.dest]:
                    args[group.title] = True
                    break

    return args


def sanity_checks(args, myid):
    """
	Check if the command line arguments are valid.

	Arguments:
	args - Command line arguments as dictionary
	myid - MPI id

	Returns:
	None
	"""

    if args["gain_file"] is not None:
        if not os.path.isfile(args["gain_file"]):
            sp_global_def.ERROR(
                "If the gain_file option is provided, the gain file must exist!",
                myid=myid,
            )

    if args["additional_dose_unadjusted"] and args["skip_dose_adjustment"]:
        sp_global_def.ERROR(
            "If the additional_dose_unadjusted option is provided, the skip_dose_adjustment cannot be provided as well!",
            myid=myid,
        )

    if args["selection_file"] is not None and not os.path.isfile(args.selection_file):
        sp_global_def.ERROR(
            "If the selection_file option is provided, the specified file needs to exist!",
            myid=myid,
        )

    if not glob.glob(args["input_micrograph_pattern"]):
        sp_global_def.ERROR(
            "No files found with the specified pattern!: {0}".format(
                args["input_micrograph_pattern"]
            ),
            myid=myid,
        )

    if os.path.exists(args["output_directory"]) and not args["overwrite"]:
        sp_global_def.ERROR("Output directory is not allowed to exist!", myid=myid)


def load_file_names_by_pattern(pattern, selection_file):
    """


	Arguments:
	args - Command line arguments as dictionary
	myid - MPI id

	Returns:
	None
	"""
    file_names_raw = sorted(glob.glob(pattern))
    sp_global_def.sxprint(
        "Found {0} movie files based on the provided pattern.".format(
            len(file_names_raw)
        )
    )

    if selection_file is not None:
        with open(selection_file, "r") as read:
            lines_raw = read.readlines()
        sp_global_def.sxprint(
            "Found {0} file_names in the selection file.".format(len(lines_raw))
        )

        selection_list = [
            os.path.basename(os.path.splitext(entry.strip())[0]) for entry in lines_raw
        ]
        file_names = [
            entry
            for entry in input_frames_raw
            if os.path.basename(os.path.splitext(entry)[0]) in selection_list
        ]
        sp_global_def.sxprint(
            "Found {0} movie files after applying the selection file.".format(
                len(file_names)
            )
        )
    else:
        file_names = file_names_raw
    return file_names


def create_unblur_command(
    input_stack,
    output_name,
    pixel_size,
    bin_factor,
    do_filter,
    kv,
    exp,
    pre_exp,
    do_expert,
    min_shift,
    outer_radius,
    b_factor,
    half_width_vert,
    half_width_hor,
    termination,
    max_iterations,
    restore_power,
    gain_file,
    first_frame,
    last_frame,
    correct_mag,
    dist_angle,
    major_scale,
    minor_scale,
):
    """Multiline Comment0"""

    # Command list
    unblur_command = []
    unblur_command.append(input_stack)
    unblur_command.append(output_name)
    unblur_command.append(pixel_size)
    unblur_command.append(bin_factor)

    if do_filter:
        unblur_command.append("yes")
        unblur_command.append(kv)
        unblur_command.append(exp)
        unblur_command.append(pre_exp)
    else:
        unblur_command.append("no")

    if do_expert:
        unblur_command.append("yes")
        unblur_command.append(min_shift)
        unblur_command.append(outer_radius)
        unblur_command.append(b_factor)
        unblur_command.append(half_width_vert)
        unblur_command.append(half_width_hor)
        unblur_command.append(termination)
        unblur_command.append(max_iterations)
        if do_filter:
            if restore_power:
                unblur_command.append("yes")
            else:
                unblur_command.append("no")
        if gain_file is None:
            unblur_command.append("yes")
        else:
            unblur_command.append("no")
            unblur_command.append(gain_file)
        unblur_command.append(first_frame)
        unblur_command.append(last_frame)

    else:
        unblur_command.append("no")

    if correct_mag:
        unblur_command.append("yes")
        unblur_command.append(dist_angle)
        unblur_command.append(major_scale)
        unblur_command.append(minor_scale)
    else:
        unblur_command.append("no")

    return "\n".join([str(entry) for entry in unblur_command])


def run():
    """
	Main function

	Arguments:
	args - Arguments as dictionary

	Returns:
	None
	"""
    args = parse_args()
    main_mpi_proc = 0
    my_mpi_proc_id = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    n_mpi_procs = mpi.mpi_comm_size(mpi.MPI_COMM_WORLD)

    # Import the file names
    sanity_checks(args, my_mpi_proc_id)
    if my_mpi_proc_id == main_mpi_proc:
        if args["Expert options"]:
            sp_global_def.sxprint(
                "Expert option detected! The program will enable expert mode!"
            )
        if args["Magnification correction"]:
            sp_global_def.sxprint(
                "Magnification correction option detected! The program will enable magnification correction mode!"
            )
        file_names = load_file_names_by_pattern(
            args["input_micrograph_pattern"], args["selection_file"]
        )
    else:
        file_names = []

    file_names = sp_utilities.wrap_mpi_bcast(file_names, main_mpi_proc)

    # Split the list indices by node
    max_proc = min(n_mpi_procs, len(file_names))
    if my_mpi_proc_id in list(range(max_proc)):
        idx_start, idx_end = sp_applications.MPI_start_end(
            len(file_names), max_proc, my_mpi_proc_id
        )
    else:
        idx_start = 0
        idx_end = 0

    nima = idx_end - idx_start
    max_nima_list = sp_utilities.wrap_mpi_gatherv(
        [nima], main_mpi_proc, mpi.MPI_COMM_WORLD
    )
    max_nima_list = sp_utilities.wrap_mpi_bcast(
        max_nima_list, main_mpi_proc, mpi.MPI_COMM_WORLD
    )
    max_nima = max(max_nima_list)
    mpi_print_id = max_nima_list.index(max_nima)

    try:
        os.makedirs(args["output_directory"])
    except OSError:
        pass
    sp_global_def.write_command(args["output_directory"])
    start_unblur = time.time()
    for idx, file_path in enumerate(file_names[idx_start:idx_end]):
        if my_mpi_proc_id == mpi_print_id:
            total_time = time.time() - start_unblur
            if idx == 0:
                average_time = 0
            else:
                average_time = old_div(total_time, float(idx))
            sp_global_def.sxprint(
                "{0: 6.2f}% => Elapsed time: {1: 6.2f}min | Estimated total time: {2: 6.2f}min | Time per micrograph: {3: 5.2f}min/mic".format(
                    old_div(100 * idx, float(max_nima)),
                    old_div(total_time, float(60)),
                    old_div((max_nima) * average_time, float(60)),
                    old_div(average_time, float(60)),
                )
            )

        file_name = os.path.basename(os.path.splitext(file_path)[0])
        file_name_out = "{0}.mrc".format(file_name)
        file_name_log = "{0}.log".format(file_name)
        file_name_err = "{0}.err".format(file_name)

        output_dir_name = os.path.join(args["output_directory"], "corrsum")
        output_dir_name_log = os.path.join(args["output_directory"], "corrsum_log")
        output_dir_name_dw = os.path.join(args["output_directory"], "corrsum_dw")
        output_dir_name_dw_log = os.path.join(
            args["output_directory"], "corrsum_dw_log"
        )
        if args["additional_dose_unadjusted"]:
            unblur_list = (
                (True, output_dir_name_dw, output_dir_name_dw_log),
                (False, output_dir_name, output_dir_name_log),
            )
        elif args["skip_dose_adjustment"]:
            unblur_list = ((False, output_dir_name, output_dir_name_log),)
        else:
            unblur_list = ((True, output_dir_name_dw, output_dir_name_dw_log),)

        for dose_adjustment, dir_name, log_dir_name in unblur_list:
            try:
                os.makedirs(dir_name)
            except OSError:
                pass
            try:
                os.makedirs(log_dir_name)
            except OSError:
                pass
            output_name = os.path.join(dir_name, file_name_out)
            output_name_log = os.path.join(log_dir_name, file_name_log)
            output_name_err = os.path.join(log_dir_name, file_name_err)
            unblur_command = create_unblur_command(
                file_path,
                output_name,
                args["pixel_size"],
                args["bin_factor"],
                dose_adjustment,
                args["voltage"],
                args["exposure_per_frame"],
                args["pre_exposure"],
                args["Expert options"],
                args["min_shift_initial"],
                args["outer_radius"],
                args["b_factor"],
                args["half_width_vert"],
                args["half_width_hor"],
                args["termination"],
                args["max_iterations"],
                bool(not args["dont_restore_noise_power"]),
                args["gain_file"],
                args["first_frame"],
                args["last_frame"],
                args["Magnification correction"],
                args["distortion_angle"],
                args["major_scale"],
                args["minor_scale"],
            )

            execute_command = r'echo "{0}" | {1}'.format(
                unblur_command, args["unblur_path"]
            )
            with open(output_name_log, "w") as log, open(output_name_err, "w") as err:
                start = time.time()
                child = subprocess.Popen(
                    execute_command, shell=True, stdout=log, stderr=err
                )
                child.wait()
                if child.returncode != 0:
                    sp_global_def.sxprint(
                        "Process failed for image {0}.\nPlease make sure that the unblur path is correct\nand check the respective logfile.".format(
                            file_path
                        )
                    )
                log.write(
                    "Time => {0:.2f} for command: {1}".format(
                        time.time() - start, execute_command
                    )
                )

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)

    if my_mpi_proc_id == mpi_print_id:
        idx = idx + 1
        total_time = time.time() - start_unblur
        average_time = old_div(total_time, float(idx))
        sp_global_def.sxprint(
            "{0: 6.2f}% => Elapsed time: {1: 6.2f}min | Estimated total time: {2: 6.2f}min | Time per micrograph: {3: 5.2f}min/mic".format(
                old_div(100 * idx, float(max_nima)),
                old_div(total_time, float(60)),
                old_div((max_nima) * average_time, float(60)),
                old_div(average_time, float(60)),
            )
        )

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.BATCH = True
    run()
    sp_global_def.BATCH = False
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()

if __name__ == "__main__":
    main()
