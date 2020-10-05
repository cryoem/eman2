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
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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
from ..libpy import sp_applications
from ..libpy import sp_global_def
from ..libpy import sp_utilities
import sys


def run():
    progname = optparse.os.path.basename(sys.argv[0])
    usage = progname + " stack_in  stack_out"
    parser = optparse.OptionParser(usage, version=sp_global_def.SPARXVERSION)
    (options, args) = parser.parse_args()

    if sp_global_def.CACHE_DISABLE:
        sp_utilities.disable_bdb_cache()

    # check length of arguments list. less than 2 is illegal
    if len(args) < 2:
        sp_global_def.sxprint("usage: " + usage)
        sp_global_def.sxprint("Please run '" + progname + " -h' for detailed options")
    # 2 is file to file copying
    elif 2 == len(args):
        # print "file to file"
        sp_applications.cpy(args[0], args[1])
    # more than 2, this means a wildcard is transformed to a list of filenams
    else:
        # print "list to file"
        #
        # XXX: note that wildcards only work for filenames! wildcards in files
        #    are processed by the shell and read in as a space-delimited list
        #    of filenames. however, this does NOT work for bdb:image1 type names,
        #    since these will not be found and processed by the shell, not being
        #    normal files. globbing of db objects will have to be done here!
        #
        # make sure we only pass the last entry of args as output file name,
        #    since [-1:] pass a list containing only the last entry....
        #
        # application.cpy
        sp_applications.cpy(args[:-1], args[-1:][0])

def main():
    sp_global_def.print_timestamp("Start")
    sp_global_def.write_command()
    run()
    sp_global_def.print_timestamp("Finish")


if __name__ == "__main__":
    main()