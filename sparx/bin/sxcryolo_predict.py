#!/usr/bin/env python
#
# Copyright (C) 2016  Thorsten Wagner (thorsten.wagner@mpi-dortmund.mpg.de)
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
import argparse
from json import dump
import subprocess

argparser = argparse.ArgumentParser(
    description='Apply crYOLO on your dataset',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

argparser.add_argument(
    'config_path',
    type=str,
    help='Specifiy the path to your config file.')

argparser.add_argument(
    'target_dir',
    type=str,
    help='Specifiy the path to your config file.')

argparser.add_argument(
    'model_path',
    type=str,
    help='Specifiy the path to your model file.')

argparser.add_argument(
    'output_dir',
    type=str,
    help='Path to where the boxfiles are written.')

argparser.add_argument(
    '--confidence_threshold',
    type=float,
    default=0.3,
    help='Confidence threshold')

argparser.add_argument(
    '--gpu',
    default=0,
    type=int,
    nargs="+",
    help="Specifiy which gpu(s) should be used. Multiple GPUs are separated by a whitespace")

argparser.add_argument(
    '--filament_mode',
    action="store_true",
    help="Specifiy if filament mode should be used")

argparser.add_argument(
    '--filament_width',
    default=100,
    type=int,
    help="Spezify the width of the filament in pixel")

argparser.add_argument(
    '--box_distance',
    default=100,
    type=int,
    help="Distance between two boxes in pixel")

argparser.add_argument(
    '--min_box_per_filament',
    default=6,
    type=int,
    help="Minimum number of boxes per filament")


def main():
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

    if type(args.gpu) is list:
        str_gpus = [str(entry) for entry in args.gpu]
    else:
        str_gpus = str(args.gpu)

    arg_gpu = ' '.join(str_gpus)

    # Run the training
    complete_command = ['cryolo_predict.py']
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
    gpu_argument = "-g=" + arg_gpu
    complete_command.append(gpu_argument)
    if do_filament_mode:
        complete_command.append("--filament")
        complete_command.append("-fw="+str(filament_width))
        complete_command.append("-mn=" + str(min_box_per_filament))
        complete_command.append("-bd=" + str(box_distance))
    subprocess.check_call(complete_command)

if __name__ == "__main__":
	main()