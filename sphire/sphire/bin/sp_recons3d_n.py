#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

import EMAN2_cppwrap
import mpi
import optparse
from ..libpy import sp_applications
from ..libpy import sp_global_def
from ..libpy import sp_utilities
import string
import sys
from builtins import range


def run():
    arglist = []
    for arg in sys.argv:
        arglist.append(arg)

    progname = optparse.os.path.basename(arglist[0])
    usage = (
        progname
        + " prj_stack volume [begin end step] --CTF --npad=ntimes_padding --list=file --group=ID --snr=SNR --sym=symmetry --verbose=(0|1) --xysize --MPI"
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)

    parser.add_option(
        "--CTF", action="store_true", default=False, help="apply CTF correction"
    )
    parser.add_option("--snr", type="float", default=1.0, help="Signal-to-Noise Ratio")
    parser.add_option("--sym", type="string", default="c1", help="symmetry")
    parser.add_option(
        "--list",
        type="string",
        help="file with list of images to be used in the first column",
    )
    parser.add_option(
        "--group",
        type="int",
        default=-1,
        help="perform reconstruction using images for a given group number (group is attribute in the header)",
    )
    parser.add_option(
        "--MPI", action="store_true", default=False, help="use MPI version "
    )
    parser.add_option(
        "--npad", type="int", default=2, help="number of times padding (default 2)"
    )
    parser.add_option(
        "--verbose",
        type="int",
        default=0,
        help="verbose level: 0 no verbose, 1 verbose",
    )
    parser.add_option(
        "--xysize", type="int", default=-1, help="user expected size at xy direction"
    )
    parser.add_option(
        "--zsize", type="int", default=-1, help="user expected size at z direction"
    )
    parser.add_option(
        "--smearstep",
        type="float",
        default=0.0,
        help="Rotational smear step (default 0.0, no smear)",
    )
    parser.add_option(
        "--interpolation_method",
        type="string",
        default="4nn",
        help="4nn, or tril: nearest neighbor, or trilinear interpolation",
    )
    parser.add_option(
        "--niter",
        type="int",
        default=10,
        help="number of iterations for iterative reconstruction",
    )
    parser.add_option(
        "--upweighted",
        action="store_true",
        default=False,
        help="apply background noise",
    )
    parser.add_option(
        "--compensate",
        action="store_true",
        default=False,
        help="compensate in reconstruction",
    )
    parser.add_option(
        "--chunk_id",
        type="int",
        default=-1,
        help="reconstruct both odd and even groups of particles",
    )
    parser.add_option(
        "--target_window_size",
        type="int",
        default=-1,
        help=" size of the targeted reconstruction ",
    )

    (options, args) = parser.parse_args(arglist[1:])

    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()

    if len(args) == 2:
        prj_stack = args[0]
        vol_stack = args[1]
        nimage = EMAN2_cppwrap.EMUtil.get_image_count(prj_stack)
        pid_list = list(range(0, nimage))
    elif len(args) == 5:
        prj_stack = args[0]
        vol_stack = args[1]
        begin = int(args[2])
        end = int(args[3])
        step = int(args[4])
        pid_list = list(range(begin, end, step))
    else:
        sp_global_def.ERROR("Incomplete list of arguments")
        return

    if options.list and options.group > -1:
        sp_global_def.ERROR("options group and list cannot be used together")
        return

    sp_global_def.BATCH = True
    if options.interpolation_method == "4nn":
        sp_applications.recons3d_n(
            prj_stack,
            pid_list,
            vol_stack,
            options.CTF,
            options.snr,
            1,
            options.npad,
            options.sym,
            options.list,
            options.group,
            options.verbose,
            options.MPI,
            options.xysize,
            options.zsize,
            options.smearstep,
            options.upweighted,
            options.compensate,
            options.chunk_id,
        )
    elif options.interpolation_method == "tril":
        if options.MPI is False:
            sp_global_def.ERROR(
                "Trilinear interpolation reconstruction has MPI version only!"
            )
            return
        sp_applications.recons3d_trl_MPI(
            prj_stack,
            pid_list,
            vol_stack,
            options.CTF,
            options.snr,
            1,
            options.npad,
            options.sym,
            options.verbose,
            options.niter,
            options.compensate,
            options.target_window_size,
        )

    else:
        sp_global_def.ERROR(
            "Wrong interpolation method. The current options are 4nn, and tril. 4nn is the default one."
        )
        return

    sp_global_def.write_command(optparse.os.path.dirname(vol_stack))
    sp_global_def.BATCH = False

def main():
    RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in optparse.os.environ
    if RUNNING_UNDER_MPI:
        mpi.mpi_init(
            0, []
        )  # On OS X, there is an error if MPI is initialized and not finalized, hence the conditional
    sp_global_def.print_timestamp("Start")
    sp_global_def.write_command()
    run()
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()
    if RUNNING_UNDER_MPI:
        mpi.mpi_finalize()

if __name__ == "__main__":
    main()
