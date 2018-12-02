#
from __future__ import print_function
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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

import EMAN2_cppwrap
import sparx_applications
import sparx_filter
import sparx_fundamentals
import sparx_global_def
import sparx_logger
import sparx_morphology
import mpi
import sparx_multi_shc
import numpy
import numpy as np
import numpy.random
import operator
import sparx_pixel_error
import sparx_projection
import scipy
import sparx_statistics
import sys
import time
import types
import sparx_utilities
pass#IMPORTIMPORTIMPORT import EMAN2
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import alignment
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import cPickle as pickle
pass#IMPORTIMPORTIMPORT import collections
pass#IMPORTIMPORTIMPORT import development
pass#IMPORTIMPORTIMPORT import filter
pass#IMPORTIMPORTIMPORT import fundamentals
pass#IMPORTIMPORTIMPORT import gc
pass#IMPORTIMPORTIMPORT import glob
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import inspect
pass#IMPORTIMPORTIMPORT import logger
pass#IMPORTIMPORTIMPORT import math
pass#IMPORTIMPORTIMPORT import morphology
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import multi_shc
pass#IMPORTIMPORTIMPORT import numpy
pass#IMPORTIMPORTIMPORT import numpy as np
pass#IMPORTIMPORTIMPORT import numpy.random
pass#IMPORTIMPORTIMPORT import operator
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import pixel_error
pass#IMPORTIMPORTIMPORT import projection
pass#IMPORTIMPORTIMPORT import random
pass#IMPORTIMPORTIMPORT import scipy
pass#IMPORTIMPORTIMPORT import socket
pass#IMPORTIMPORTIMPORT import sparx
pass#IMPORTIMPORTIMPORT import statistics
pass#IMPORTIMPORTIMPORT import subprocess
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import time
pass#IMPORTIMPORTIMPORT import types
pass#IMPORTIMPORTIMPORT import utilities
from builtins import range
pass#IMPORTIMPORTIMPORT from global_def import *
pass#IMPORTIMPORTIMPORT import numpy.random


#  06-12-14 code lifted
'''
def Numrinit(first_ring, last_ring, skip=1, mode="F"):
	#  This is to test equal length rings
	"""This function calculates the necessary information for the 2D 
	   polar interpolation. For each ring, three elements are recorded:
	   numr[i*3]:  Radius of this ring
	   numr[i*3+1]: Total number of samples of all inner rings+1
	   		(Or, the beginning point of this ring)
	   numr[i*3+2]: Number of samples of this ring. This number is an 
	   		FFT-friendly power of the 2.
			
	   "F" means a full circle interpolation
	   "H" means a half circle interpolation
	"""
	MAXFFT = 32768
	pass#IMPORTIMPORTIMPORT from math import pi

	if (mode == 'f' or mode == 'F'): dpi = 2*pi
	else:                            dpi = pi
	numr = []
	lcirc = 1
	#  This is for testing equal length rings
	ip = 128
	for k in xrange(first_ring, last_ring+1, skip):
		numr.append(k)
		numr.append(lcirc)
		numr.append(ip)
		lcirc += ip		
	return  numr
"""Multiline Comment6"""
			# The following code is used when mirror is not considered
			retvals = Util.Crosrng_e(crefim, cimage, numr, 0, 0.0)
			qn = retvals["qn"]
			if qn >= peak:
				sx = -ix
				sy = -iy
				ang = ang_n(retvals["tot"], mode, numr[-1])
				peak = qn
				mirror = 0
			"""Multiline Comment7"""
