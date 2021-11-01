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
# Author: Pawel A.Penczek (Pawel.A.Penczek@uth.tmc.edu)
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

import optparse
from sphire.libpy import sp_applications
from sphire.libpy import sp_global_def
from sphire.libpy import sp_utilities
import sys


def run():
    arglist = []
    for arg in sys.argv:
        arglist.append(arg)

    progname = optparse.os.path.basename(arglist[0])
    usage = (
        progname
        + " stack --params='parm1 parm2 parm3 ...' --zero --one --set=number --randomize --rand_alpha --import=file --export=file --print --backup --suffix --restore --delete"
    )
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)

    parser.add_option("--params", type="string", default=None, help="parameter list")
    parser.add_option(
        "--zero", action="store_true", default=False, help="set parameter to zero"
    )
    parser.add_option(
        "--one", action="store_true", default=False, help="set parameter to one"
    )
    parser.add_option(
        "--set",
        type="float",
        default=0.0,
        help="set parameter to a value (different from 0.0)",
    )
    parser.add_option(
        "--randomize",
        action="store_true",
        default=False,
        help="set parameter to randomized value",
    )
    parser.add_option(
        "--rand_alpha",
        action="store_true",
        default=False,
        help="set all angles to randomized value",
    )
    parser.add_option(
        "--import",
        type="string",
        dest="fimport",
        default=None,
        help="import parameters from file",
    )
    parser.add_option(
        "--export",
        type="string",
        dest="fexport",
        default=None,
        help="export parameters to file",
    )
    parser.add_option(
        "--print",
        action="store_true",
        dest="fprint",
        default=False,
        help="print parameters",
    )
    parser.add_option(
        "--backup", action="store_true", default=False, help="backup parameters"
    )
    parser.add_option(
        "--suffix",
        type="string",
        default="_backup",
        help="suffix for xform name in backup",
    )
    parser.add_option(
        "--restore", action="store_true", default=False, help="restore parameters"
    )
    parser.add_option(
        "--delete", action="store_true", default=False, help="delete parameters"
    )
    parser.add_option(
        "--consecutive",
        action="store_true",
        default=False,
        help="set selected parameter to consecutive integers starting from 0",
    )
    parser.add_option(
        "--list",
        type="string",
        default=None,
        help="Indices list containing the same amount of rows as the import file",
    )

    (options, args) = parser.parse_args(arglist[1:])

    if not options.fprint:
        sp_global_def.print_timestamp("Start")
        sp_global_def.write_command()

    if len(args) != 1:
        sp_global_def.sxprint("Usage: " + usage)
        sp_global_def.ERROR(
            "Invalid number of parameters provided. Please see usage information above."
        )
        return

    if options.params == None:
        sp_global_def.sxprint("Usage: " + usage)
        sp_global_def.ERROR(
            "No parameters provided. Please see usage information above."
        )
        return

    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()
    sp_applications.header(
        args[0],
        options.params,
        options.zero,
        options.one,
        options.set,
        options.randomize,
        options.rand_alpha,
        options.fimport,
        options.fexport,
        options.fprint,
        options.backup,
        options.suffix,
        options.restore,
        options.delete,
        options.consecutive,
        options.list,
    )
    if not options.fprint:
        sp_global_def.print_timestamp("Finish")

def main():
    run()

if __name__ == "__main__":
    main()
