#!/usr/bin/env python
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
import mpi
import optparse
from ..libpy import sp_global_def
from ..libpy import sp_logger
from ..libpy import sp_multi_shc
from ..libpy import sp_user_functions
from ..libpy import sp_utilities
import sys
from builtins import range





def run(args):

    progname = optparse.os.path.basename(sys.argv[0])
    usage = (
        progname
        + " stack  [output_directory] --ir=inner_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1 --ref_a=S --sym=c1"
    )
    usage += """

stack			2D images in a stack file: (default required string)
directory		output directory name: into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results) (default required string)
"""

    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    parser.add_option(
        "--radius",
        type="int",
        default=29,
        help="radius of the particle: has to be less than < int(nx/2)-1 (default 29)",
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
        "--mask3D", type="string", default=None, help="3D mask file: (default sphere)"
    )
    parser.add_option(
        "--moon_elimination",
        type="string",
        default="",
        help="elimination of disconnected pieces: two arguments: mass in KDa and pixel size in px/A separated by comma, no space (default none)",
    )
    parser.add_option(
        "--ir",
        type="int",
        default=1,
        help="inner radius for rotational search: > 0 (default 1)",
    )

    # 'radius' and 'ou' are the same as per Pawel's request; 'ou' is hidden from the user
    # the 'ou' variable is not changed to 'radius' in the 'sparx' program. This change is at interface level only for sxviper.
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option("--ou", type="int", default=-1, help=optparse.SUPPRESS_HELP)

    parser.add_option(
        "--rs",
        type="int",
        default=1,
        help="step between rings in rotational search: >0 (default 1)",
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
        "--nruns",
        type="int",
        default=6,
        help="GA population: aka number of quasi-independent volumes (default 6)",
    )
    parser.add_option(
        "--doga",
        type="float",
        default=0.1,
        help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga: (default 0.1)",
    )
    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--npad",
        type="int",
        default=2,
        help="padding size for 3D reconstruction (default=2)",
    )
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
        "--debug",
        action="store_true",
        default=False,
        help="debug info printout: (default False)",
    )

    ##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    parser.add_option(
        "--return_options",
        action="store_true",
        dest="return_options",
        default=False,
        help=optparse.SUPPRESS_HELP,
    )

    parser.add_option("--filament_width", default=-1)
    # parser.add_option("--an",       type="string", default= "-1",               help="NOT USED angular neighborhood for local searches (phi and theta)")
    # parser.add_option("--CTF",      action="store_true", default=False,         help="NOT USED Consider CTF correction during the alignment ")
    # parser.add_option("--snr",      type="float",  default= 1.0,                help="NOT USED Signal-to-Noise Ratio of the data (default 1.0)")
    # (options, args) = parser.parse_args(sys.argv[1:])

    required_option_list = ["radius"]
    (options, args) = parser.parse_args(args)
    # option_dict = vars(options)
    # print parser

    if options.return_options:
        return parser

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
            sp_global_def.ERROR("Missing parameter. Please see above")
            return

    if len(args) < 2 or len(args) > 3:
        sp_global_def.sxprint("Usage: " + usage)
        sp_global_def.sxprint("Please run '" + progname + " -h' for detailed options")
        sp_global_def.ERROR(
            "Invalid number of parameters used. Please see usage information above."
        )
        return

    log = sp_logger.Logger(sp_logger.BaseLogger_Files())

    # 'radius' and 'ou' are the same as per Pawel's request; 'ou' is hidden from the user
    # the 'ou' variable is not changed to 'radius' in the 'sparx' program. This change is at interface level only for sxviper.
    options.ou = options.radius
    runs_count = options.nruns
    mpi_rank = mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD)
    mpi_size = mpi.mpi_comm_size(
        mpi.MPI_COMM_WORLD
    )  # Total number of processes, passed by --np option.

    if mpi_rank == 0:
        all_projs = EMAN2_cppwrap.EMData.read_images(args[0])
        subset = list(range(len(all_projs)))
        # if mpi_size > len(all_projs):
        # 	ERROR('Number of processes supplied by --np needs to be less than or equal to %d (total number of images) ' % len(all_projs), 'sxviper', 1)
        # 	mpi.mpi_finalize()
        # 	return
    else:
        all_projs = None
        subset = None

    outdir = args[1]
    error = 0
    if mpi_rank == 0:
        if mpi_size % options.nruns != 0:
            sp_global_def.ERROR(
                "Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.",
                action=0,
            )
            error = 1

        if optparse.os.path.exists(outdir):
            sp_global_def.ERROR(
                "Output directory '%s' exists, please change the name and restart the program"
                % outdir,
                action=0,
            )
            error = 1
        sp_global_def.LOGFILE = optparse.os.path.join(outdir, sp_global_def.LOGFILE)

    mpi.mpi_barrier(mpi.MPI_COMM_WORLD)
    error = sp_utilities.bcast_number_to_all(
        error, source_node=0, mpi_comm=mpi.MPI_COMM_WORLD
    )
    if error == 1:
        return

    if mpi_rank == 0:
        optparse.os.makedirs(outdir)
        sp_global_def.write_command(outdir)

    if outdir[-1] != "/":
        outdir += "/"
    log.prefix = outdir

    # if len(args) > 2:
    # 	ref_vol = get_im(args[2])
    # else:
    # ref_vol = None

    options.user_func = sp_user_functions.factory[options.function]

    options.CTF = False
    options.snr = 1.0
    options.an = -1.0
    options.filament_width = -1
    out_params, out_vol, out_peaks = sp_multi_shc.multi_shc(
        all_projs, subset, runs_count, options, mpi_comm=mpi.MPI_COMM_WORLD, log=log
    )

def main():
    mpi.mpi_init(0, [])
    sp_global_def.print_timestamp("Start")
    main(sys.argv[1:])
    sp_global_def.print_timestamp("Finish")
    mpi.mpi_finalize()

if __name__ == "__main__":
    main()