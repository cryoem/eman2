#!/usr/bin/env python

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


from __future__ import print_function
from __future__ import division
import argparse
from ..libpy import sp_global_def
import subprocess


argparser = argparse.ArgumentParser(
    description="Run automatic 2d class selelection (Cinderella)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)


argparser.add_argument("cinderella_path", default=None, type=str, help="Path to your ")

argparser.add_argument("inputstack", help="Path to your class stack (.mrcs/.hdf)")

argparser.add_argument(
    "output_dir", type=str, help="Path to where the boxfiles are written."
)

argparser.add_argument(
    "model_path", type=str, help="Specifiy the path to your model file."
)


argparser.add_argument(
    "-t",
    "--confidence_threshold",
    default=0.5,
    type=float,
    help="Classes with a confidence higher as that threshold are classified as good.",
)

argparser.add_argument("--gpu", default=-1, type=int, help="GPU to run on.")

argparser.add_argument(
    "-b",
    "--batch_size",
    default=32,
    type=int,
    help="Number of mini-batches during prediction.",
)


def run():
    args = argparser.parse_args()
    cindy_path = args.cinderella_path
    input_stack = args.inputstack
    out_dir = args.output_dir
    model_path = args.model_path
    conf_thresh = args.confidence_threshold
    batch_size = args.batch_size
    gpu = args.gpu

    complete_command = [cindy_path]

    arg_input_stack = "-i=" + str(input_stack)
    complete_command.append(arg_input_stack)

    arg_out_dir = "-o=" + str(out_dir)
    complete_command.append(arg_out_dir)

    arg_model_path = "-w=" + str(model_path)
    complete_command.append(arg_model_path)

    arg_conf_thresh = "-t=" + str(conf_thresh)
    complete_command.append(arg_conf_thresh)

    arg_batch_size = "-b=" + str(batch_size)
    complete_command.append(arg_batch_size)

    if gpu != -1:
        arg_gpu = "--gpu=" + str(gpu)
        complete_command.append(arg_gpu)

    sp_global_def.write_command()
    subprocess.check_call(complete_command)
    sp_global_def.write_command(out_dir)

def main():
    sp_global_def.print_timestamp("Start")
    run()
    sp_global_def.print_timestamp("Finish")

if __name__ == "__main__":
    main()
