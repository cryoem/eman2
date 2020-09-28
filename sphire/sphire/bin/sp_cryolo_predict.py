#!/usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2016-2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
#
# Copyright (c) 2016-2019 Max Planck Institute of Molecular Physiology
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

from __future__ import print_function
from __future__ import division
import argparse
from ..libpy import sp_global_def
import subprocess


argparser = argparse.ArgumentParser(
    description="Apply crYOLO on your dataset",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

argparser.add_argument(
    "config_path", type=str, help="Specifiy the path to your config file."
)

argparser.add_argument(
    "target_dir", type=str, help="Specifiy the path to your config file."
)

argparser.add_argument(
    "model_path", type=str, help="Specifiy the path to your model file."
)

argparser.add_argument(
    "output_dir", type=str, help="Path to where the boxfiles are written."
)

argparser.add_argument("--cryolo_predict_path", default=None, type=str)

argparser.add_argument(
    "--confidence_threshold", type=float, default=0.3, help="Confidence threshold"
)

argparser.add_argument(
    "--gpu",
    default="0",
    type=str,
    help="Specifiy which gpu(s) should be used. Multiple GPUs are separated by a comma",
)

argparser.add_argument("--otf", action="store_true", help="On-the-fly filtering")

argparser.add_argument(
    "--min_distance",
    default=0,
    type=int,
    help="Particles with a distance less than this value (in pixel) will be removed",
)

argparser.add_argument(
    "--filament_mode",
    action="store_true",
    help="Specifiy if filament mode should be used",
)

argparser.add_argument(
    "--filament_width",
    default=100,
    type=int,
    help="Spezify the width of the filament in pixel",
)

argparser.add_argument(
    "--box_distance",
    default=-1,
    type=int,
    help="Distance between two filament boxes in pixel",
)

argparser.add_argument(
    "--min_box_per_filament",
    default=6,
    type=int,
    help="Minimum number of boxes per filament",
)

argparser.add_argument(
    "--nosplit",
    action="store_true",
    help="(FILAMENT MODE) The filament mode does not split to curved filaments",
)

argparser.add_argument(
    "--nomerging",
    action="store_true",
    help="(FILAMENT MODE) The filament mode does not merge filaments",
)

argparser.add_argument(
    "--gpu_fraction",
    type=float,
    default=1.0,
    help="Specify the fraction of memory per GPU used by crYOLO during prediction. Only values between 0.0 and 1.0 are allowed.",
)

argparser.add_argument(
    "-nc",
    "--num_cpu",
    type=int,
    default=-1,
    help="Number of CPUs used during prediction. By default it will use half of the available CPUs.",
)


def run():
    # Read arguments
    args = argparser.parse_args()

    config_path = args.config_path
    target_dir = args.target_dir
    model_path = args.model_path
    output_dir = args.output_dir
    confidence_threshold = args.confidence_threshold

    do_filament_mode = args.filament_mode
    filament_width = args.filament_width
    min_box_per_filament = args.min_box_per_filament
    box_distance = args.box_distance
    gpu_fraction = 1.0
    if args.gpu_fraction < 1.0 and args.gpu_fraction > 0.0:
        gpu_fraction = args.gpu_fraction
    num_cpu = args.num_cpu

    str_gpus = [str(entry).strip() for entry in args.gpu.split(",")]
    arg_gpu = " ".join(str_gpus)

    no_merging = args.nomerging
    no_split = args.nosplit
    cryolo_predict_path = args.cryolo_predict_path

    # Run the training
    complete_command = ["cryolo_predict.py"]
    if cryolo_predict_path is not None:
        complete_command = [cryolo_predict_path]

    config_argument = "-c=" + str(config_path)
    complete_command.append(config_argument)
    weights_argument = "-w=" + str(model_path)
    complete_command.append(weights_argument)
    input_argument = "-i=" + str(target_dir)
    complete_command.append(input_argument)
    output_argument = "-o=" + str(output_dir)
    complete_command.append(output_argument)
    thresh_argument = "-t=" + str(confidence_threshold)
    complete_command.append(thresh_argument)

    if args.min_distance > 0:
        min_dist_arg = "-d=" + str(int(args.min_distance))
        complete_command.append(min_dist_arg)

    gpu_argument = "-g " + arg_gpu
    complete_command.append(gpu_argument)
    gpu_fraction_arg = "--gpu_fraction=" + str(gpu_fraction)
    complete_command.append(gpu_fraction_arg)
    if num_cpu != -1:
        num_cpu_arg = "--num_cpu=" + str(num_cpu)
        complete_command.append(num_cpu_arg)

    if args.otf:
        complete_command.append("--otf")

    if do_filament_mode:
        complete_command.append("--filament")
        complete_command.append("-fw=" + str(filament_width))
        complete_command.append("-mn=" + str(min_box_per_filament))
        if box_distance > 0:
            complete_command.append("-bd=" + str(box_distance))
        if no_merging:
            complete_command.append("--nomerging")
        if no_split:
            complete_command.append("--nosplit")
    sp_global_def.write_command()
    subprocess.check_call(complete_command)
    sp_global_def.write_command(output_dir)


def main():
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")

if __name__ == "__main__":
    main()
