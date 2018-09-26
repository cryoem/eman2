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

def main():
    # Read arguments
    args = argparser.parse_args()

    config_path = args.config_path
    target_dir = args.target_dir
    model_path = args.model_path
    output_dir = args.output_dir
    confidence_threshold = args.confidence_threshold

    if type(args.gpu) is list:
        str_gpus = [str(entry) for entry in args.gpu]
    else:
        str_gpus = str(args.gpu)

    arg_gpu = ' '.join(str_gpus)

    # Run the training
    config_argument = "-c=" + str(config_path)
    weights_argument = "-w=" + str(model_path)
    input_argument = "-i=" + str(target_dir)
    output_argument = "-o=" + str(output_dir)
    thresh_argument = "-t=" + str(confidence_threshold)
    gpu_argument = "-g=" + arg_gpu
    subprocess.check_call(['cryolo_predict.py', config_argument, weights_argument, input_argument, output_argument,thresh_argument,gpu_argument])

if __name__ == "__main__":
	main()