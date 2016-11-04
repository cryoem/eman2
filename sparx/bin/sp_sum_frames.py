#! /usr/bin/env python
#
# Copyright (C) 2016  Markus Stabrin (markus.stabrin@mpi-dortmund.mpg.de)
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

from sparx import EMData, EMUtil
from sys import argv
from global_def import SPARXVERSION
from optparse import OptionParser, SUPPRESS_HELP
import global_def

parser = OptionParser('', version=SPARXVERSION)
(options, args) = parser.parse_args(argv[1:])

global_def.BATCH = True

# Get the size of the input files
hdr = EMData(args[0])
nx = hdr['nx']
ny = hdr['ny']
nz = hdr['nz']
print(nx, ny, nz)
hdr2 = hdr.get_clip(nz)


outvol = EMData(nx, ny, nz)


