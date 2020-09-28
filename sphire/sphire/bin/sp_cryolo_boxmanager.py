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
    description="crYOLO boxmanager",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

argparser.add_argument(
    "--cryolo_bm_path", type=str, help="Specifiy cryolo boxmanager executeable."
)

argparser.add_argument(
    "--target_dir", type=str, help="Specifiy the path to your image directory."
)

argparser.add_argument(
    "--box_dir", type=str, help="Specifiy the path to a box file directory."
)


def run():
    # Read arguments
    args = argparser.parse_args()

    target_dir = args.target_dir
    box_dir = args.box_dir

    call = ["cryolo_boxmanager.py"]
    if args.cryolo_bm_path:
        call = [args.cryolo_bm_path]
    if target_dir:
        input_argument = "-i=" + str(target_dir)
        call.append(input_argument)
        if box_dir:
            box_argument = "-b=" + str(box_dir)
            call.append(box_argument)
    print(call)
    subprocess.check_call(call)

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.write_command()
    run()
    sp_global_def.print_timestamp("Finish")
if __name__ == "__main__":
    main()
