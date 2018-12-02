#
from __future__ import print_function
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holfds
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

import EMAN2_cppwrap
import sparx_alignment
import sparx_applications
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_isac
import mpi
import numpy
import os
import sparx_pixel_error
import random
import sparx_statistics
import time
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import isac
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pixel_error
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import sparx
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
pass#IMPORTIMPORTIMPORT from global_def import SPARX_MPI_TAG_UNIVERSAL

def generate_random_averages(data, K, rand_seed = -1):
	pass#IMPORTIMPORTIMPORT from random import shuffle, seed, randint
	#  I prefer to take random images....  PAP
	if rand_seed == -1:  random.seed(random.randint(1,2000111222))
	else:                random.seed(rand_seed)
	ndata = len(data)
	ll = list(range(ndata))
	random.shuffle(ll)
	return [data[ll[i]].copy() for i in range(K)]

	"""Multiline Comment5"""
	return avgs



